"""
Module to Determine share of shared coast between eez and administrative units.

Based on `scripts/shared_coast.py` from  https://github.com/calliope-project/solar-and-wind-potentials
"""

from functools import partial
from multiprocessing import Pool

import geopandas as gpd
import pandas as pd
from eurocalliopelib.geo import EPSG3035


def allocate_eezs(
    path_to_eez: str,
    path_to_units: str,
    path_to_continental_units: str,
    threads: int,
    polygon_area_share_threshold: float,
    resolution: str,
    path_to_output: str,
):
    """
    Allocate Exclusive Economic Zones (EEZs) to Euro-Calliope units.
    If "continental" or "national" resolution, the EEZs will be applied directly to units.
    If any subnational resolution exists, the national EEZs will be shared between those
    units which have a coast, based on the length of coast.
    This is a rough estimate and is impacted by the number of small islands just off the
    coast (which increase a unit's share of the coast artificially). We mitigate this
    effect slightly by ignoring all islands below a certain size (using `polygon_area_share_threshold`).

    Args:
        path_to_eez (str):
            Path to shapefile containing Exclusive Economic Zones. Must include the
            columns "iso_ter1" (country IDs), "mrgid" (EEZ IDs), and "geometry".
        path_to_units (str):
            Path to shapefile containing Euro-Calliope units at specfied "resolution".
        path_to_continental_units (str):
            Path to shapefile containing Euro-Calliope units at the continental
            resolution.
        threads (int):
            Maximum number of threads allowed for parallelisation of operations.
        polygon_area_share_threshold (float):
            Fraction of a Euro-Calliope unit's area below which a component of the unit
            can be ignored (e.g. small islands).
        resolution (str):
            Resolution of the Euro-Calliope units in "path_to_units".
        path_to_output (str):
            Path to save pandas.Series with fraction of coast of each Euro-Calliope unit
            ("id") for every EEZ polygon ("mrgid")
    """
    eez = gpd.read_file(path_to_eez).to_crs(EPSG3035)
    units = gpd.read_file(path_to_units).to_crs(EPSG3035)

    # ASSUME: we can ignore EEZs shared by two countries
    # (these are very small far offshore, so removing them does not impact the results)
    eez = eez[eez.pol_type.str.lower() != "joint regime"]

    if resolution == "continental":
        share = pd.Series(
            index=pd.MultiIndex.from_product(
                (["EUR"], eez.mrgid.unique()), names=["id", "mrgid"]
            ),
            data=1,
        )

    elif resolution == "national":
        share = pd.Series(
            index=eez.set_index(["iso_ter1", "mrgid"]).index.rename(["id", "mrgid"]),
            data=1,
        )
    else:
        eez.geometry = eez.geometry.map(_buffer_if_necessary)
        continental_units = gpd.read_file(path_to_continental_units).to_crs(EPSG3035)
        with Pool(int(threads)) as pool:
            share_of_coast_length = pool.map(
                partial(
                    _share_of_coast_length,
                    eez=eez,
                    units=units,
                    continental_units=continental_units,
                    polygon_area_share_threshold=polygon_area_share_threshold,
                ),
                eez.iso_ter1.unique(),
            )
        share = pd.concat(share_of_coast_length)

    share.index = share.index.set_levels(
        levels=(
            share.index.levels[0].map(lambda x: x.replace(".", "-")),
            share.index.levels[1].astype(int),
        ),
        level=[0, 1],
    )

    share_sum = share.groupby(level="mrgid").sum().reindex(eez.mrgid.unique()).fillna(0)
    assert (
        (share_sum > 0.99) & (share_sum < 1.01)  # ensure we haven't missed any area
    ).all(), share_sum

    share.rename("shared_coast_fraction").to_csv(path_to_output)


def _get_coastal_units_as_linestrings(
    units: gpd.GeoDataFrame,
    continental_units: gpd.GeoDataFrame,
    polygon_area_share_threshold: float,
):
    """
    Get the outline of all Euro-Calliope units which sit on the coast
    (i.e., will have some share of the EEZ assigned to them)
    """
    # slightly increase the sub-continental polygon size so that those on the coast
    # slightly overlap the continent polygon boundary.
    units["geometry"] = units.geometry.buffer(units.total_bounds.mean() * 1e-6)
    non_coastal_units = gpd.sjoin(units, continental_units, op="within")
    coastal_units = units.loc[~units.id.isin(non_coastal_units.id_left)].set_index("id")

    # Simplify geometries to get rid of tiny islands that slow down the computation and
    # artificially increase one unit's share of the total length of coast.
    simplified_coastal_geometries = _simplify_geometries(
        coastal_units, polygon_area_share_threshold
    )

    # Return the boundary of the polygons
    return gpd.GeoDataFrame(
        geometry=simplified_coastal_geometries.boundary, crs=EPSG3035
    )


def _simplify_geometries(units, polygon_area_share_threshold):
    """
    Remove tiny islands from units to speed up the later intersection.
    Any polygons in a multipolygon A with an area below
    (polygon_area_share_threshold * A) will be removed.
    """
    all_polygons = units.geometry.explode()
    return (
        all_polygons.where(
            all_polygons.area.div(all_polygons.area.groupby(level="id").sum())
            > polygon_area_share_threshold
        )
        .dropna()
        .reset_index()
        .dissolve("id")
        .reindex(units.index)
    )


def _share_of_coast_length(
    country, eez, units, continental_units, polygon_area_share_threshold
):
    """
    Parallelisable sub-function which allocates a share of EEZ units ("mrgid") connected
    to a specfic country ("iso_ter1") to Euro-Calliope units ("id") in that same country
    ("country_code").
    """
    coastal_unit_boundaries = _get_coastal_units_as_linestrings(
        units[units.country_code == country].copy(),
        continental_units,
        polygon_area_share_threshold,
    )
    unit_intersection = gpd.overlay(
        coastal_unit_boundaries.reset_index(),
        eez.loc[eez.iso_ter1 == country, ["mrgid", "geometry"]],
        how="intersection",
    )
    coast_length_ratio = (
        unit_intersection.set_index(["id", "mrgid"])
        .length.groupby("mrgid")
        .transform(lambda x: x / x.sum())
    )

    return coast_length_ratio


# TODO: replace with shapely.make_valid (requires updating many geo.yaml env dependencies to udpate to shapely 1.8.2)
def _buffer_if_necessary(shape):
    """Fix the basins shapes which are invalid.

    Following the advice given here:
    https://github.com/Toblerity/Shapely/issues/344
    """
    if not shape.is_valid:
        shape = shape.buffer(0.0)
    assert shape.is_valid
    return shape


if __name__ == "__main__":
    allocate_eezs(
        path_to_eez=snakemake.input.eez,
        path_to_units=snakemake.input.units,
        path_to_continental_units=snakemake.input.continental_units,
        threads=int(snakemake.threads),
        polygon_area_share_threshold=snakemake.params.polygon_area_share_threshold,
        resolution=snakemake.wildcards.resolution,
        path_to_output=snakemake.output[0],
    )
