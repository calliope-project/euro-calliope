import math

import numpy as np
import pandas as pd
import geopandas as gpd
from rasterio.transform import from_origin
from rasterstats import zonal_stats
from shapely.geometry import Point


ID_DTYPE = np.uint16
ID_NO_DATA_VALUE = 64001
INDEX_EPSILON = 10e-3
DEPRECATED_GRID_SIZE_IN_M = 50000 # old style capacity factors are on a grid of 50km size

WGS84_EPSG = 4326
WGS84 = f"EPSG:{WGS84_EPSG}"
EPSG3035 = "EPSG:3035"


def area_weighted_time_series(shapes, spatiotemporal):
    """Forms area weighted time series for collections of shapes.

    Inputs:
        * locations: a GeoDataFrame of shapes, each will receive one time series
        * spatiotemporal: a DataArray with dimensions "x", "y", "timestep" in CRS WGS84
    """
    assert_correct_form(shapes, spatiotemporal)
    id_map, ts = transform_to_int_indexed(spatiotemporal)
    weighted_ts_ids_per_shape = zonal_stats(
        shapes.geometry,
        id_map,
        affine=infer_transform(spatiotemporal),
        categorical=True,
        all_touched=True,
        nodata=ID_NO_DATA_VALUE
    )

    return pd.DataFrame(
        index=ts.timestep.to_index(),
        data={
            shape_id: shape_time_series(weighted_ts_ids, ts)
            for shape_id, weighted_ts_ids in zip(shapes.index, weighted_ts_ids_per_shape)
        }
    )


def assert_correct_form(shapes, spatiotemporal):
    assert shapes.crs
    assert "crs" in spatiotemporal.attrs.keys()
    assert shapes.crs == gpd.tools.crs.CRS(spatiotemporal.attrs["crs"])
    assert "y" in spatiotemporal.dims, "Expect dimension 'y'"
    assert "x" in spatiotemporal.dims, "Expect dimension 'x'"
    assert "timestep" in spatiotemporal.dims, "Expect dimension 'timestep'"


def shape_time_series(weighted_ts_ids, ts):
    ts_ids = [int(idx) for idx in weighted_ts_ids.keys()]
    relevant_ts = (
        ts
        .sel({"z": ts_ids})
        .copy()
        .dropna("z", how="any") # remove points without data
    )
    weights = (
        pd
        .Series(index=ts_ids, data=weighted_ts_ids.values())
        .reindex(relevant_ts.z)
        .transform(lambda x: x / x.sum())
        .rename_axis(index="z")
        .to_xarray()
    )
    weighted_cfs = relevant_ts * weights
    return weighted_cfs.sum(dim="z").to_series()


def transform_to_int_indexed(spatiotemporal):
    """Create a raster map of ids to timeseries.

    Each point on the map links to a timeseries.
    """
    x_min = spatiotemporal.x.min().item()
    x_max = spatiotemporal.x.max().item()
    y_min = spatiotemporal.y.min().item()
    y_max = spatiotemporal.y.max().item()
    resolution = infer_resolution(spatiotemporal)
    width = (x_max - x_min) / resolution + 1
    height = (y_max - y_min) / resolution + 1
    assert isclose(round(width), width) # diff is purely numerics
    assert isclose(round(height), height) # diff is purely numerics
    width = round(width)
    height = round(height)
    raster = np.ones(shape=(height, width), dtype=ID_DTYPE) * ID_NO_DATA_VALUE
    stacked_spatiotemporal = spatiotemporal.stack(z=["x", "y"])
    for n, z in enumerate(stacked_spatiotemporal.z):
        x, y = z.item()
        index_x = (x - x_min) / resolution
        index_y = (y_max - y) / resolution
        assert isclose(round(index_x), index_x) # diff is purely numerics
        assert isclose(round(index_y), index_y) # diff is purely numerics
        int_index_x = round(index_x)
        int_index_y = round(index_y)
        raster[int_index_y, int_index_x] = n
    stacked_spatiotemporal = (
        stacked_spatiotemporal
        .assign_coords(z=range(len(stacked_spatiotemporal.z))) # convert to int index
    )
    return raster, stacked_spatiotemporal


def isclose(a, b):
    return math.isclose(a, b, abs_tol=INDEX_EPSILON, rel_tol=0)


def infer_resolution(spatiotemporal):
    y_diffs = spatiotemporal.y.diff("y")
    x_diffs = spatiotemporal.x.diff("x")
    resolution_x = x_diffs[0].item()
    assert (x_diffs == resolution_x).all()
    resolution_y = y_diffs[0].item()
    assert (y_diffs == resolution_y).all()
    assert resolution_x == resolution_y
    return resolution_x


def infer_transform(spatiotemporal):
    resolution = infer_resolution(spatiotemporal)
    x_min = spatiotemporal.x.min()
    y_max = spatiotemporal.y.max()
    return from_origin(
        west=x_min - resolution / 2,
        north=y_max + resolution / 2,
        xsize=resolution,
        ysize=resolution
    )


def convert_old_style_capacity_factor_time_series(ts):
    """DEPRECATED: Converts published capacity factor data to new format.

    Existing capacity factor time series that are used in this workflow are published here:
    https://zenodo.org/record/3899687

    The old format is integer-indexed so that the integer refers to the location. `lat` and `lon`
    are variables of this dataset (they should be part of the index, but are not). In addition,
    while the data have been created by forming an equally-spaced grid on a projected plane
    (EPSG:3035), the data themselves are given in WGS84 and thus not on an equally-spaced grid.
    We will not use this format anymore in the future. Thus, we need to support it only as long
    as we did not update the published data.

    The new format is coordinate-indexed so that `timestep`, `x`, and `y` are the dimensions of
    the data. The disadvantage of this format is that it potentially includes a lot of `nans` in
    locations in which no data exist. As both netcdf4 and xarray support sparse datasets, this
    should not be a big disadvantage.

    This function takes data in the old format and converts them to the new format. The function is
    deprecated and should be removed as soon as the published data is updated with the new format.
    """
    site_id_map = {
        site_id.item(): (ts["lon"].sel(site_id=site_id).item(), ts["lat"].sel(site_id=site_id).item())
        for site_id in ts.site_id
    }
    gdf = (
        gpd
        .GeoDataFrame(
            data={"site_id": [site_id for site_id in site_id_map.keys()]},
            geometry=[Point(lon, lat) for lon, lat in site_id_map.values()],
            crs=WGS84
        )
        .set_index("site_id")
        .to_crs(EPSG3035)
    )
    gdf["x"] = [round(point.coords[0][0], ndigits=0) for point in gdf.geometry] # round to meter
    gdf["y"] = [round(point.coords[0][1], ndigits=0) for point in gdf.geometry] # round to meter
    ts["x"] = gdf["x"]
    ts["y"] = gdf["y"]
    ts = ts.set_index(site_id=["x", "y"]).unstack("site_id")
    ts = ts["electricity"].rename(time="timestep")

    # make sure the resolution is uniform
    x = ts.x
    y = ts.y
    assert ((x.diff("x") == DEPRECATED_GRID_SIZE_IN_M).sum() / x.count()).item() > 0.9 # most are fine
    assert ((y.diff("y") == DEPRECATED_GRID_SIZE_IN_M).sum() / x.count()).item() > 0.9 # most are fine
    ts = ts.reindex(
        x=np.arange(x[0], x[-1] + 1, step=DEPRECATED_GRID_SIZE_IN_M),
        y=np.arange(y[0], y[-1] + 1, step=DEPRECATED_GRID_SIZE_IN_M)
    )
    ts.attrs["crs"] = EPSG3035
    return ts
