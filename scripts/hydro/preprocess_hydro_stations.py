import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import pycountry

from eurocalliopelib import utils


MAXIMUM_NUMBER_OF_DROPPED_STATIONS = 5


def preprocess_stations(path_to_stations, country_codes, path_to_basins, buffer_size,
                        path_to_phs_storage_capacities, scale_phs, path_to_output):
    stations = pd.read_csv(path_to_stations, index_col=0)
    stations = drop_kosovo(stations)
    stations = drop_stations_without_installed_capacity(stations)
    stations = fix_station_country_code(stations)
    stations = resolve_index_duplicates(stations)
    stations = fill_missing_storage_capacity_values(stations)
    if scale_phs:
        stations = scale_phs_acording_to_geth(stations, path_to_phs_storage_capacities)
    stations = drop_out_of_scope(stations, country_codes)
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


def drop_out_of_scope(stations, country_codes):
    return stations.where(stations.country_code.isin(country_codes)).dropna(axis="index", subset=["country_code"])


def fill_missing_storage_capacity_values(stations):
    for mask in [stations.type == "HDAM", stations.type == "HPHS"]:
        stations.loc[mask, "storage_capacity_MWh"] = estimate_missing_storage_capacity_values(stations.loc[mask])
    return stations


def estimate_missing_storage_capacity_values(stations):
    # ASSUME country specific median E/P ratio for missing values, global median where no country specific available
    e_to_p = stations.storage_capacity_MWh / stations.installed_capacity_MW
    based_on_global = e_to_p.median() * stations.installed_capacity_MW
    based_on_country_specific = stations.merge(
        e_to_p.groupby(stations.country_code).median().rename("country_specific"),
        left_on="country_code",
        right_index=True,
        how="left"
    ).loc[:, "country_specific"] * stations.installed_capacity_MW
    storage_capacity_MWh = (
        stations.storage_capacity_MWh
        .where(pd.notnull, other=based_on_country_specific)
        .where(pd.notnull, other=based_on_global)
    )
    assert not storage_capacity_MWh.isnull().any()
    return storage_capacity_MWh


def scale_phs_acording_to_geth(stations, path_to_geth_national_capacities):
    """Scale PHS storage capacities to match (Geth et al., 2015).

    Storage capacities of pumped hydro within the JRC database may seem too high.
    Thus, they are scaled here so that the national numbers of (Geth et al., 2015) are fulfilled.

    Geth, F., Brijs, T., Kathan, J., Driesen, J., Belmans, R., 2015. An overview of large-scale
    stationary electricity storage plants in Europe: Current status and new developments. Renewable
    and Sustainable Energy Reviews 52, 1212–1227. https://doi.org/10.1016/j.rser.2015.07.145
    """
    hphs_mask = stations.type == "HPHS"
    storage_capacity_share = (
        stations
        .loc[hphs_mask]
        .groupby("country_code")
        .storage_capacity_MWh
        .transform(lambda x: x / x.sum())
        .rename("storage_capacity_share")
    )
    geth_national_capacities = (
        read_national_phs_storage_capacities(path_to_geth_national_capacities)
        .reindex(stations.loc[hphs_mask].country_code)
    )
    geth_national_capacities.index = stations.loc[hphs_mask].index
    new_storage_capacities = storage_capacity_share * geth_national_capacities
    stations.loc[hphs_mask, "storage_capacity_MWh"] = new_storage_capacities
    return stations


def read_national_phs_storage_capacities(path_to_data):
    data = pd.read_csv(path_to_data, index_col=0)
    data.index = data.index.map(utils.eu_country_code_to_iso3)
    return (
        data["storage-capacity-gwh"]
        .mul(1000)
        .rename("storage_capacity_mwh")
    )


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
    ill_placed_stations = stations[~is_in_basin]
    for station_id, station in ill_placed_stations.iterrows():
        new_point = new_coords(station, buffer_size, hydrobasins)
        stations.loc[station_id, "lon"] = new_point.coords[0][0]
        stations.loc[station_id, "lat"] = new_point.coords[0][1]
    assert stations.apply(station_in_any_basin(hydrobasins), axis="columns").all()
    return stations


def new_coords(station, buffer_size, hydrobasins):
    point = Point(station.lon, station.lat)
    basin_id = hydrobasins.distance(point).idxmin()
    closest_basin = hydrobasins.loc[basin_id]
    return point.buffer(buffer_size).intersection(closest_basin.geometry).representative_point()


if __name__ == "__main__":
    preprocess_stations(
        path_to_stations=snakemake.input.stations,
        country_codes=[pycountry.countries.lookup(country).alpha_3
                       for country in snakemake.params.countries],
        path_to_phs_storage_capacities=snakemake.input.phs_storage_capacities,
        scale_phs=snakemake.params.scale_phs,
        path_to_basins=snakemake.input.basins,
        buffer_size=snakemake.params.buffer_size,
        path_to_output=snakemake.output[0]
    )
