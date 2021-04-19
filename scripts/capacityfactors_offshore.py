"""Generate offshore capacityfactor time series."""
import pandas as pd
import geopandas as gpd
import xarray as xr

from eurocalliopelib.geo import area_weighted_time_series, convert_old_style_capacity_factor_time_series

EPSG3035 = "EPSG:3035"


def capacityfactors(path_to_eez, path_to_shared_coast, path_to_timeseries,
                    threshold, area_weighting_resolution_in_m, path_to_result, year=None):
    """Generate offshore capacityfactor time series for each location."""
    eez = gpd.read_file(path_to_eez).set_index("id").to_crs(EPSG3035).geometry
    shared_coast = pd.read_csv(path_to_shared_coast, index_col=0)
    shared_coast.index = shared_coast.index.map(lambda x: x.replace(".", "-"))
    ts = xr.open_dataset(path_to_timeseries)
    if year:
        ts = ts.sel(time=str(year))
    ts = convert_old_style_capacity_factor_time_series(ts)

    capacityfactors_per_eez = area_weighted_time_series(
        shapes=eez,
        spatiotemporal=ts,
        resolution=area_weighting_resolution_in_m
    )
    capacityfactors = _allocate_to_onshore_locations(capacityfactors_per_eez, shared_coast)
    capacityfactors.where(capacityfactors >= threshold, 0).to_csv(path_to_result)


def _allocate_to_onshore_locations(capacityfactors_per_eez, shared_coast):
    return pd.DataFrame(
        index=capacityfactors_per_eez.index,
        data={
            location_id: _onshore_timeseries(location_id, capacityfactors_per_eez, shared_coast)
            for location_id in shared_coast.index
        }
    )


def _onshore_timeseries(location_id, capacityfactors_per_eez, shared_coast):
    weights = shared_coast.loc[location_id, :].transform(lambda x: x / x.sum())
    return (capacityfactors_per_eez * weights).sum(axis="columns")


if __name__ == '__main__':
    capacityfactors(
        path_to_eez=snakemake.input.eez,
        path_to_shared_coast=snakemake.input.shared_coast,
        path_to_timeseries=snakemake.input.timeseries,
        threshold=float(snakemake.params.threshold),
        year=snakemake.params.year if snakemake.params.trim_ts else None,
        area_weighting_resolution_in_m=snakemake.params.area_weighting_resolution_in_m,
        path_to_result=snakemake.output[0]
    )
