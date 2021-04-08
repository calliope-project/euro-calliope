import importlib

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

MAXIMUM_NUMBER_OF_DROPPED_STATIONS = 5
# import utility module
assert "utils" in snakemake.input.keys(), "Make sure utils.py module is in your snakemake inputs."
spec = importlib.util.spec_from_file_location("utils", snakemake.input.utils)
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


def preprocess_stations(path_to_stations, path_to_basins, buffer_size, path_to_output):
    stations = pd.read_csv(path_to_stations, index_col=0)
    stations = resolve_index_duplicates(stations)
    stations = drop_kosovo(stations)
    stations = drop_stations_without_installed_capacity(stations)
    stations = fix_station_country_code(stations)
    stations = ensure_stations_within_basins(stations, path_to_basins, buffer_size)
    stations.to_csv(
        path_to_output,
        header=True,
        index=True
    )


def resolve_index_duplicates(stations):
    stations.index = stations.index.where(~stations.index.duplicated(), stations.index + '_b')
    assert not stations.index.duplicated().any()
    return stations.rename_axis(index="id")


def drop_kosovo(stations):
    return stations.where(stations.country_code != "XK").dropna(axis="index", subset=["country_code"])


def drop_stations_without_installed_capacity(stations):
    assert stations.installed_capacity_MW.isna().sum() <= MAXIMUM_NUMBER_OF_DROPPED_STATIONS
    return stations.dropna(subset=["installed_capacity_MW"])


def fix_station_country_code(stations):
    stations["country_code"] = stations.country_code.apply(utils.eu_country_code_to_iso3)
    return stations


def station_in_any_basin(basins):
    def station_in_any_basin(station):
        point = Point(station.lon, station.lat)
        return basins.geometry.intersects(point).sum() > 0
    return station_in_any_basin


def ensure_stations_within_basins(stations, path_to_basins, buffer_size):
    hydrobasins = gpd.read_file(path_to_basins)
    is_in_basin = stations.apply(station_in_any_basin(hydrobasins), axis="columns")
    stations = stations.copy()
    ill_placed_stations = stations[~is_in_basin]
    for station_id, station in ill_placed_stations.iterrows():
        new_point = new_coords(station, buffer_size, hydrobasins)
        stations.loc[station_id, "lon"] = new_point.coords[0][0]
        stations.loc[station_id, "lat"] = new_point.coords[0][1]
    assert stations[~is_in_basin].apply(station_in_any_basin(hydrobasins), axis="columns").all()
    return stations


def new_coords(station, buffer_size, hydrobasins):
    point = Point(station.lon, station.lat)
    basin_id = hydrobasins.distance(point).idxmin()
    closest_basin = hydrobasins.loc[basin_id]
    return point.buffer(buffer_size).intersection(closest_basin.geometry).representative_point()


if __name__ == "__main__":
    preprocess_stations(
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        buffer_size=snakemake.params.buffer_size,
        path_to_output=snakemake.output[0]
    )
