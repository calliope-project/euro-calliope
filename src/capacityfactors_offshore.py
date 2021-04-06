"""Generate offshore capacityfactor time series."""
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import xarray as xr

SITE_ID_DIM = "site_id"
WEIGHT_VAR = "weight"
CAPACITY_FACTOR_VAR = "electricity"


def capacityfactors(path_to_eez, path_to_shared_coast, path_to_id_map, path_to_timeseries,
                    threshold, path_to_result, year=None):
    """Generate offshore capacityfactor time series for each location."""
    eez = gpd.read_file(path_to_eez).set_index("id").geometry
    shared_coast = pd.read_csv(path_to_shared_coast, index_col=0)
    shared_coast.index = shared_coast.index.map(lambda x: x.replace(".", "-"))
    ts = xr.open_dataset(path_to_timeseries)
    if year:
        ts = ts.sel(time=str(year))
    with rasterio.open(path_to_id_map, "r") as src:
        id_map = src.read(1)
        transform = src.transform
        nodata = src.nodata

    capacityfactors_per_eez = _capacityfactors(eez, id_map, transform, nodata, ts)
    capacityfactors = _allocate_to_onshore_locations(capacityfactors_per_eez, shared_coast)
    capacityfactors.where(capacityfactors >= threshold, 0).to_csv(path_to_result)


def _capacityfactors(locations, id_map, transform, nodata, ts):
    weighted_ts_ids_per_location = zonal_stats(
        locations.geometry,
        id_map,
        affine=transform,
        categorical=True,
        all_touched=True,
        nodata=nodata
    )
    return pd.DataFrame(
        index=ts.time.to_index(),
        data={
            location_id: _location_time_series(weighted_ts_ids, ts)
            for location_id, weighted_ts_ids in zip(locations.index, weighted_ts_ids_per_location)
        }
    )


def _location_time_series(weighted_ts_ids, ts):
    ts_ids = list(weighted_ts_ids.keys())
    ts_weights = pd.Series(list(weighted_ts_ids.values())).transform(lambda x: x / x.sum()).values
    relevant_ts = ts.sel({SITE_ID_DIM: ts_ids}).copy()
    relevant_ts[WEIGHT_VAR] = ((SITE_ID_DIM), ts_weights)
    weighted_cfs = relevant_ts[CAPACITY_FACTOR_VAR] * relevant_ts[WEIGHT_VAR]
    return weighted_cfs.sum(dim=SITE_ID_DIM).to_series()


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
        path_to_id_map=snakemake.input.ids,
        path_to_timeseries=snakemake.input.timeseries,
        threshold=float(snakemake.params.threshold),
        year=snakemake.params.year if snakemake.params.trim_ts else None,
        path_to_result=snakemake.output[0]
    )
