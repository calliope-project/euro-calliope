import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import pycountry

MAXIMUM_NUMBER_OF_DROPPED_STATIONS = 5


def preprocess_stations(path_to_stations, path_to_basins, buffer_size, path_to_output):
    stations = pd.read_csv(path_to_stations, index_col=0)
    stations = resolve_index_duplicates(stations)
    stations = drop_kosovo(stations)
    stations = drop_stations_without_installed_capacity(stations)
    stations = fix_station_country_code(stations)
    hydrobasins = gpd.read_file(path_to_basins)
    stations = ensure_stations_in_basins(stations, hydrobasins, buffer_size)
    stations.to_csv(
        path_to_output,
        header=True,
        index=True
    )


def resolve_index_duplicates(stations):
    stations = stations.copy()
    new_index = stations.index.values
    mask = stations.index.duplicated()
    new_index[mask] = [idx + "_b" for idx in stations.index[mask]]
    stations.index = new_index
    assert not stations.index.duplicated().any()
    return stations.rename_axis(index="id")


def drop_kosovo(stations):
    return stations.where(stations.country_code != "XK").dropna(axis="index", subset=["country_code"])


def drop_stations_without_installed_capacity(stations):
    assert stations.installed_capacity_MW.isna().sum() <= MAXIMUM_NUMBER_OF_DROPPED_STATIONS
    return stations.dropna(subset=["installed_capacity_MW"])


def fix_station_country_code(stations):
    stations["country_code"] = stations.country_code.apply(eu_country_code_to_iso3)
    return stations


def station_in_any_basin(basins):
    def station_in_any_basin(station):
        point = Point(station.lon, station.lat)
        return basins.geometry.intersects(point).sum() > 0
    return station_in_any_basin


def ensure_stations_in_basins(stations, hydrobasins, buffer_size):
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


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(eu_country_code) == 2, "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)
    if eu_country_code.lower() == "el":
        iso2 = "gr"
    elif eu_country_code.lower() == "uk":
        iso2 = "gb"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3


if __name__ == "__main__":
    preprocess_stations(
        path_to_stations=snakemake.input.stations,
        path_to_basins=snakemake.input.basins,
        buffer_size=snakemake.params.buffer_size,
        path_to_output=snakemake.output[0]
    )
