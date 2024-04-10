import math

import geopandas as gpd
import pandas as pd
import rasterio
import xarray as xr
from rasterstats import zonal_stats

EPSG_3035 = "EPSG:3035"
WGS84 = "EPSG:4326"
GRIDBOX_SIZE = 25000  # MERRA-2 grid
idx = pd.IndexSlice


# FIXME this will most likely be entirely replaced once a more generalised way to deal
# with spatial units is available
def gridded_weather_population(
    path_to_population,
    path_to_locations,
    path_to_wind_speed,  # FIXME: this is used to get the coordinates, but nothing else
    lat_name,
    lon_name,
    out_path,
):
    # preprocessed weather data is gridded across Europe with WGS84 projection and
    # includes hourly temperature (in K), average wind speed (in m/s) at every grid point.
    wind_speed_ds = xr.open_dataset(path_to_wind_speed)

    # Locations are shapefiles at the resolution of interest (e.g. countries for resolution `national`)
    locations = gpd.read_file(path_to_locations)

    gridbox_points = gpd.GeoDataFrame(
        geometry=gpd.points_from_xy(
            wind_speed_ds[lon_name].values, wind_speed_ds[lat_name].values
        ),
        index=wind_speed_ds.site.to_index(),
        crs=WGS84,
    )
    # To go from a grid of points to a grid of boxes filling the entire space,
    # we `buffer` to create a grid of circles whose edges just touch,
    # then we define the envelope of that circle to create a square.
    gridbox_points = gridbox_points.to_crs(EPSG_3035)
    gridbox = gpd.GeoDataFrame(
        gridbox_points.index.to_frame(),
        geometry=gridbox_points.buffer(GRIDBOX_SIZE).envelope,
    )
    # `Overlay` creates new shapes that are either complete gridboxes or partial ones that
    # sit inside a specific location.
    gridboxes_mapped_to_locations = gpd.overlay(
        gridbox.to_crs(WGS84), locations.to_crs(WGS84)
    )

    with rasterio.open(path_to_population) as src:
        array = src.read(1)
        meta = src.meta
        affine = src.transform

        population_per_complete_or_partial_gridbox_polygon = zonal_stats(
            gridboxes_mapped_to_locations.to_crs(meta["crs"]),
            array,
            affine=affine,
            stats="sum",
            nodata=meta["nodata"],
        )
        gridboxes_mapped_to_locations["population"] = [
            i["sum"] for i in population_per_complete_or_partial_gridbox_polygon
        ]

        population_per_country = zonal_stats(
            locations.to_crs(meta["crs"]),
            array,
            affine=affine,
            stats="sum",
            nodata=meta["nodata"],
        )
        locations["population"] = [i["sum"] for i in population_per_country]

    # Confirm that the total population is valid (i.e. we haven't picked up or lost regions).
    # This is a test that the gridboxes cover all land with population that we are interested in.
    assert math.isclose(
        locations.population.sum(),
        gridboxes_mapped_to_locations.population.sum(),
        abs_tol=10**3,
    )

    population_da = xr.DataArray.from_series(
        gridboxes_mapped_to_locations.set_index("site").population
    )
    population_da.to_netcdf(out_path)


if __name__ == "__main__":
    gridded_weather_population(
        path_to_population=snakemake.input.population,
        path_to_locations=snakemake.input.locations,
        path_to_wind_speed=snakemake.input.wind_speed,
        lat_name=snakemake.params.lat_name,
        lon_name=snakemake.params.lon_name,
        out_path=snakemake.output[0],
    )
