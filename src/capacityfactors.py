"""Generate capacityfactor time series."""
import pandas as pd
import geopandas as gpd
import rasterio
from rasterstats import zonal_stats
import xarray as xr

SITE_ID_DIM = "site_id"
WEIGHT_VAR = "weight"
CAPACITY_FACTOR_VAR = "electricity"


def capacityfactors(path_to_locations, path_to_id_map, path_to_timeseries, path_to_result):
    """Generate capacityfactor time series for each location."""
    locations = gpd.read_file(path_to_locations).set_index("id").geometry
    ts = xr.open_dataset(path_to_timeseries)
    with rasterio.open(path_to_id_map, "r") as src:
        id_map = src.read(1)
        transform = src.transform
        nodata = src.nodata

    capacityfactors = _capacityfactors(locations, id_map, transform, nodata, ts)
    capacityfactors.to_csv(path_to_result)


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
        index=ts.time,
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


if __name__ == '__main__':
    capacityfactors(
        path_to_locations=snakemake.input.locations,
        path_to_id_map=snakemake.input.ids,
        path_to_timeseries=snakemake.input.timeseries,
        path_to_result=snakemake.output[0]
    )
