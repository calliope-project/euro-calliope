"""Generate offshore capacityfactor time series."""

import geopandas as gpd
import pandas as pd
import xarray as xr
from eurocalliopelib.geo import (
    area_weighted_time_series,
    convert_old_style_capacity_factor_time_series,
)

EPSG3035 = "EPSG:3035"


def capacityfactors(
    path_to_eez,
    path_to_shared_coast,
    path_to_timeseries,
    cf_threshold,
    path_to_result,
    gridcell_overlap_threshold,
    first_year=None,
    final_year=None,
):
    """Generate offshore capacityfactor time series for each location."""
    eez = gpd.read_file(path_to_eez).set_index("mrgid").to_crs(EPSG3035).geometry
    shared_coast = pd.read_csv(path_to_shared_coast, index_col=0)
    shared_coast.index = shared_coast.index.map(lambda x: x.replace(".", "-"))
    shared_coast.columns = shared_coast.columns.astype(int)

    ts = xr.open_dataset(path_to_timeseries)
    ts = ts.sel(time=slice(first_year, final_year))
    # xarray will silently miss the fact that data doesn't exist with slice
    if (
        first_year is not None
        and first_year not in ts.time.to_index().year.astype(str).unique()
    ):
        raise ValueError(
            f"Cannot access capacity factor data for timeseries {path_to_timeseries} "
            f"with a start year of {first_year}."
        )
    if (
        final_year is not None
        and final_year not in ts.time.to_index().year.astype(str).unique()
    ):
        raise ValueError(
            f"Cannot access capacity factor data for timeseries {path_to_timeseries} "
            f"with an end year of {final_year}."
        )
    ts = convert_old_style_capacity_factor_time_series(ts)

    capacityfactors_per_eez = area_weighted_time_series(
        shapes=eez,
        spatiotemporal=ts,
        gridcell_overlap_threshold=gridcell_overlap_threshold,
    )
    capacityfactors = _allocate_to_onshore_locations(
        capacityfactors_per_eez, shared_coast
    )
    capacityfactors.where(capacityfactors >= cf_threshold, 0).to_csv(path_to_result)


def _allocate_to_onshore_locations(capacityfactors_per_eez, shared_coast):
    return pd.DataFrame(
        index=capacityfactors_per_eez.index,
        data={
            location_id: _onshore_timeseries(
                location_id, capacityfactors_per_eez, shared_coast
            )
            for location_id in shared_coast.index
        },
    )


def _onshore_timeseries(location_id, capacityfactors_per_eez, shared_coast):
    weights = shared_coast.loc[location_id, :].transform(lambda x: x / x.sum())
    return (capacityfactors_per_eez * weights).sum(axis="columns")


if __name__ == "__main__":
    capacityfactors(
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_timeseries=snakemake.input.timeseries,
        cf_threshold=float(snakemake.params.cf_threshold),
        gridcell_overlap_threshold=float(snakemake.params.gridcell_overlap_threshold),
        first_year=str(snakemake.params.first_year)
        if snakemake.params.trim_ts
        else None,
        final_year=str(snakemake.params.final_year)
        if snakemake.params.trim_ts
        else None,
        path_to_result=snakemake.output[0],
    )
