from pathlib import Path

import atlite
import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point
import pycountry


def determine_water_inflow(path_to_cutout, path_to_stations, path_to_basins, path_to_output):
    path_to_cutout = Path(path_to_cutout)
    plants = read_plants(path_to_stations)

    inflow_m3 = water_inflow(plants, path_to_cutout, path_to_basins)
    (xr.merge([plants.to_xarray(), inflow_m3])
       .drop("geometry")
       .to_netcdf(path_to_output))


def read_plants(path_to_stations):
    plants = pd.read_csv(path_to_stations, index_col=0)
    plants["country_code"] = plants["country_code"].map(lambda iso2: pycountry.countries.lookup(iso2).alpha_3)
    plants = plants[plants["type"].isin(["HROR", "HDAM"])]
    return gpd.GeoDataFrame(
        plants,
        geometry=list(map(Point, zip(plants.lon, plants.lat)))
    )


def water_inflow(plants, path_to_cutout, path_to_basins):
    cutout = atlite.Cutout(
        name=path_to_cutout.name,
        cutout_dir=path_to_cutout.parent
    )
    inflow = (cutout.hydro(plants, path_to_basins)
                    .ffill(dim="time") # FIXME probably not at all necessary
                    .bfill(dim="time") # FIXME only necessary because 2015 data not yet downloaded
                    .rename(plant="id")
                    .rename("inflow_m3"))
    # TODO fix negative values earlier
    # Inflow data can contain negative values, because runoff data contains tiny negative values.
    # see https://github.com/FRESNA/atlite/issues/17
    # see https://gist.github.com/timtroendle/832ea97fc3f594560ce7dcaf95af4fe5
    # It would be better to fix those issues in the runoff rather than the inflow data.
    return inflow.where(inflow >= 0, 0)


if __name__ == "__main__":
    determine_water_inflow(
        path_to_cutout=snakemake.input.runoff,
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        path_to_output=snakemake.output[0]
    )
