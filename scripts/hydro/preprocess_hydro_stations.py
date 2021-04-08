import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

MAXIMUM_NUMBER_OF_DROPPED_STATIONS = 5


def preprocess_stations(path_to_stations, path_to_basins, buffer_size, path_to_output):
    stations = pd.read_csv(path_to_stations, index_col=0)
    assert stations.installed_capacity_MW.isna().sum() <= MAXIMUM_NUMBER_OF_DROPPED_STATIONS
    stations.dropna(subset=["installed_capacity_MW"], inplace=True)
    hydrobasins = gpd.read_file(path_to_basins)
    is_in_basin = stations.apply(station_in_any_basin(hydrobasins), axis="columns")
    stations = move_stations_into_basins(stations, hydrobasins, is_in_basin, buffer_size)
    assert stations[~is_in_basin].apply(station_in_any_basin(hydrobasins), axis="columns").all()
    stations.to_csv(
        path_to_output,
        header=True,
        index=True
    )


def station_in_any_basin(basins):
    def station_in_any_basin(station):
        point = Point(station.lon, station.lat)
        return basins.geometry.intersects(point).sum() > 0
    return station_in_any_basin


def move_stations_into_basins(stations, hydrobasins, is_in_basin, buffer_size):
    stations = stations.copy()
    ill_placed_stations = stations[~is_in_basin]
    for station_id, station in ill_placed_stations.iterrows():
        new_point = new_coords(station, buffer_size, hydrobasins)
        stations.loc[station_id, "lon"] = new_point.coords[0][0]
        stations.loc[station_id, "lat"] = new_point.coords[0][1]
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
