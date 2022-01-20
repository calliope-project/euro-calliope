import atlite
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point
import pycountry


def determine_water_inflow(paths_to_cutout, path_to_stations, path_to_basins, first_year, final_year, path_to_output):
    plants = read_plants(path_to_stations)

    inflow_m3 = water_inflow(plants, paths_to_cutout, path_to_basins)
    (xr.merge([plants.to_xarray(), inflow_m3])
       .drop("geometry")
       .sel(time=slice(str(first_year), str(final_year)))
       .to_netcdf(path_to_output))


def read_plants(path_to_stations):
    plants = pd.read_csv(path_to_stations, index_col=0)
    plants["country_code"] = plants["country_code"].map(lambda iso2: pycountry.countries.lookup(iso2).alpha_3)
    plants = plants[plants["type"].isin(["HROR", "HDAM"])]
    return gpd.GeoDataFrame(
        plants,
        geometry=list(map(Point, zip(plants.lon, plants.lat)))
    )


def water_inflow(plants, paths_to_cutout, path_to_basins):
    cutout = atlite.Cutout(data=xr.open_mfdataset(paths_to_cutout, concat_dim="time"))
    inflow = (cutout.hydro(plants, path_to_basins)
                    .rename(plant="id")
                    .rename("inflow_m3"))
    return inflow


if __name__ == "__main__":
    determine_water_inflow(
        path_to_cutout=snakemake.input.runoff,
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        path_to_output=snakemake.output[0]
    )
