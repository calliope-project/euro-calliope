import numpy as np
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point
from geopandas.tools import overlay


DEPRECATED_GRID_SIZE_IN_M = 50000 # old style capacity factors are on a grid of 50km size

WGS84_EPSG = 4326
WGS84 = f"EPSG:{WGS84_EPSG}"
EPSG3035 = "EPSG:3035"


def area_weighted_time_series(shapes, spatiotemporal, gridcell_overlap_threshold):
    """Forms area weighted time series for collections of shapes.

    Inputs:
        * locations: a GeoDataFrame of shapes, each will receive one time series
        * spatiotemporal: a DataArray with dimensions "x", "y", "timestep"
    """
    assert_correct_form(shapes, spatiotemporal)
    shapes = shapes.rename_axis(index="shape_id").reset_index()
    stacked_spatiotemporal = spatiotemporal.stack(xy=["x", "y"])
    weights_and_values = xr.Dataset({
        "value": stacked_spatiotemporal,
        "weight": weights_between_shape_and_xy(shapes, stacked_spatiotemporal)
    })
    return pd.DataFrame(
        index=spatiotemporal.timestep.to_index(),
        data={
            shape_id: weighted_time_series(
                weights_and_values.sel(shape_id=shape_id), gridcell_overlap_threshold
            )
            for shape_id in shapes.shape_id
        }
    )


def assert_correct_form(shapes, spatiotemporal):
    assert shapes.crs
    assert "crs" in spatiotemporal.attrs.keys()
    assert shapes.crs == gpd.tools.crs.CRS(spatiotemporal.attrs["crs"])
    assert "y" in spatiotemporal.dims, "Expect dimension 'y'"
    assert "x" in spatiotemporal.dims, "Expect dimension 'x'"
    assert "timestep" in spatiotemporal.dims, "Expect dimension 'timestep'"


def weights_between_shape_and_xy(shapes, stacked_spatiotemporal):
    x, y = zip(*stacked_spatiotemporal.xy.values)
    grid_gdf = (
        gpd
        .GeoSeries(
            gpd.points_from_xy(x=x, y=y, crs=stacked_spatiotemporal.crs),
            crs=shapes.crs
        )
        .buffer(infer_resolution(stacked_spatiotemporal.unstack("xy")) / 2)
        .envelope
    )
    index_gdf = gpd.GeoDataFrame(
        geometry=grid_gdf,
        data={"xy": stacked_spatiotemporal.xy},
        crs=shapes.crs
    )
    overlaid = overlay(index_gdf, shapes, how="intersection")
    overlaid["area"] = overlaid.area
    overlaid["weight"] = overlaid.groupby("shape_id").area.transform(lambda area: area / area.sum())
    weights = (
        overlaid
        .set_index(["xy", "shape_id"])
        .loc[:, "weight"]
        .to_xarray()
        .reindex_like(stacked_spatiotemporal) # add all xy's even without overlay
        .fillna(0) # xy's without overlay have 0 weight
    )
    return weights


def weighted_time_series(weights_and_values, gridcell_overlap_threshold):
    ds = (  # drop all locations with weight == 0 or value == np.nan
        weights_and_values
        .where(weights_and_values.weight > 0)
        .dropna(subset=["weight", "value"], dim="xy", how="all")
    )
    assert ds.weight.sum() >= gridcell_overlap_threshold
    if ds.weight.sum().round(5) != 1:
        print(
            f"Weight of shape_id {ds.shape_id.item()} only adds up to {ds.weight.sum().item()}. "
            "Scaling to a sum of 1."
        )
        ds["weight"] = ds.weight / ds.weight.sum()
    return (ds * ds.weight).value.sum("xy", skipna=False)


def infer_resolution(spatiotemporal):
    y_diffs = abs(spatiotemporal.y.diff("y"))
    x_diffs = abs(spatiotemporal.x.diff("x"))
    resolution_x = x_diffs[0].item()
    assert (x_diffs == resolution_x).all()
    resolution_y = y_diffs[0].item()
    assert (y_diffs == resolution_y).all()
    assert resolution_x == resolution_y, "Resolutions in x and y must be equal."
    return resolution_x


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
