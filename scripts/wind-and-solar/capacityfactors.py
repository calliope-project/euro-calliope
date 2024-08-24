"""Generate capacityfactor time series."""

import geopandas as gpd
import xarray as xr
from eurocalliopelib.geo import (
    EPSG3035,
    area_weighted_time_series,
    convert_old_style_capacity_factor_time_series,
)


def capacityfactors(
    path_to_locations,
    path_to_timeseries,
    path_to_timeseries_with_coordinates,
    cf_threshold,
    path_to_result,
    gridcell_overlap_threshold,
    first_year=None,
    final_year=None,
):
    """Generate capacityfactor time series for each location."""
    locations = (
        gpd.read_file(path_to_locations).set_index("id").to_crs(EPSG3035).geometry
    )
    locations.index = locations.index.map(lambda x: x.replace(".", "_"))
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
    if ("lat" not in ts.dims) or (
        "lon" not in ts.dims
    ):  # rooftop pv comes without coordinates
        ts_with_latlon = xr.open_dataset(path_to_timeseries_with_coordinates)
        ts["lat"] = ts_with_latlon["lat"]
        ts["lon"] = ts_with_latlon["lon"]
    ts = convert_old_style_capacity_factor_time_series(ts)

    capacityfactors = area_weighted_time_series(
        shapes=locations,
        spatiotemporal=ts,
        gridcell_overlap_threshold=gridcell_overlap_threshold,
    )
    capacityfactors.where(capacityfactors >= cf_threshold, 0).rename_axis(
        index="timesteps"
    ).to_csv(path_to_result)


if __name__ == "__main__":
    capacityfactors(
        path_to_locations=snakemake.input.locations,
        path_to_timeseries=snakemake.input.timeseries,
        path_to_timeseries_with_coordinates=snakemake.input.coordinates,
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
