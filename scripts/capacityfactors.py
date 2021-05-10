"""Generate capacityfactor time series."""
import geopandas as gpd
import xarray as xr

from eurocalliopelib.geo import area_weighted_time_series, convert_old_style_capacity_factor_time_series, EPSG3035


def capacityfactors(path_to_locations, path_to_timeseries, path_to_timeseries_with_coordinates, threshold,
                    path_to_result, year=None):
    """Generate capacityfactor time series for each location."""
    locations = gpd.read_file(path_to_locations).set_index("id").to_crs(EPSG3035).geometry
    locations.index = locations.index.map(lambda x: x.replace(".", "-"))
    ts = xr.open_dataset(path_to_timeseries)
    if year:
        ts = ts.sel(time=str(year))
    if ("lat" not in ts.dims) or ("lon" not in ts.dims): # rooftop pv comes without coordinates
        ts_with_latlon = xr.open_dataset(path_to_timeseries_with_coordinates)
        ts["lat"] = ts_with_latlon["lat"]
        ts["lon"] = ts_with_latlon["lon"]
    ts = convert_old_style_capacity_factor_time_series(ts)

    capacityfactors = area_weighted_time_series(
        shapes=locations,
        spatiotemporal=ts
    )
    capacityfactors.where(capacityfactors >= threshold, 0).to_csv(path_to_result)


if __name__ == '__main__':
    capacityfactors(
        path_to_locations=snakemake.input.locations,
        path_to_timeseries=snakemake.input.timeseries,
        path_to_timeseries_with_coordinates=snakemake.input.coordinates,
        threshold=float(snakemake.params.threshold),
        year=snakemake.params.year if snakemake.params.trim_ts else None,
        path_to_result=snakemake.output[0]
    )
