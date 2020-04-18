import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

INVALID_ENTRIES = [
    "H2328" # ASSUME "Parteen Weir": the numbers don't match, and it seems to be a duplicate of Ardnacrucha (H585)
]
ROMANIAN_PHS = pd.DataFrame( # Capacity according to (Geth et al., 2015), all in the same region.
    index=["HXXX"],
    data={
        "name": "Geth 2015 PHS",
        "installed_capacity_MW": 200,
        "type": "HPHS",
        "country_code": "RO",
        "lat": 44.052607,
        "lon": 24.584081,
        "storage_capacity_MWh": 10.2
    }
).rename_axis(index="id")


def filter_stations(path_to_stations, path_to_basins, path_to_output):
    stations = pd.read_csv(path_to_stations, index_col=0)
    stations = add_romanian_phs(stations)
    hydrobasins = gpd.read_file(path_to_basins)
    is_in_basin = stations.apply(station_in_any_basin(hydrobasins), axis="columns")
    stations[is_in_basin].drop(index=INVALID_ENTRIES).to_csv(
        path_to_output,
        header=True,
        index=True
    )


def add_romanian_phs(stations):
    return pd.concat(
        [stations, ROMANIAN_PHS],
        axis=0,
        sort=False
    )


def station_in_any_basin(basins):
    def station_in_any_basin(station):
        point = Point(station.lon, station.lat)
        return basins.geometry.intersects(point).sum() > 0
    return station_in_any_basin


if __name__ == "__main__":
    filter_stations(
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        path_to_output=snakemake.output[0]
    )
