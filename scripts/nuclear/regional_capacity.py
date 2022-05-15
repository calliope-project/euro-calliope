import pandas as pd
import shapely.geometry
import geopandas as gpd
import pycountry


def regionalise_nuclear_capacity(
    path_to_power_plant_database: str, path_to_units: str,
    nuclear_capacity_scenario: str, countries: list, resolution: str,
    path_to_output: str
):
    """
    Use current geolocations of nuclear capacity in Europe to assign nuclear
    capacity to Euro-Calliope regions.

    Args:
        path_to_power_plant_database (str):
            location of JRC powerplant database "OPEN_UNITS" CSV file.
        path_to_units (str):
            location of euro-calliope regions geojson/shapefile.
        nuclear_capacity_scenario (str):
            nuclear scenario name or (if scenario != 'current') path to scenario file with
            minimum and maximum capacities of future nuclear capacity per country in MW
            (must include column headers [country, min, max]).
        countries (list):
            List of valid country names to filter any future nuclear capacity data on.
            Country names should match official full names that can be looked up by the
            pycountry package (e.g. `France`, `United Kingdom`).
        resolution (str):
            Euro-Calliope spatial resolution.
        path_to_output (str):
            location of output CSV file containing nuclear capacities
            (either exact - 'equals' - or 'min' and 'max') per Euro-Calliope region.
    """
    power_plants = pd.read_csv(path_to_power_plant_database)
    units = gpd.read_file(path_to_units)

    capacity_current = power_plants[power_plants.type_g == 'Nuclear']
    capacity_current_points = [
        shapely.geometry.Point(xy)
        for xy in zip(capacity_current.lon, capacity_current.lat)
    ]
    capacity_current_gdf = gpd.GeoDataFrame(
        capacity_current, geometry=capacity_current_points, crs='EPSG:4326'
    )
    capacity_current_per_region = (
        gpd.overlay(capacity_current_gdf.to_crs(units.crs), units)
        .groupby(["id", "country_code"]).sum()
        .capacity_g  # Generating unit capacity, net (MW)
    )

    if nuclear_capacity_scenario == "current":
        capacity_per_region = (
            capacity_current_per_region
            .to_frame("installed_capacity_nuclear_equals_MW")
        )
    else:
        future_capacity = _get_future_capacity_from_config_file(
            nuclear_capacity_scenario, resolution, countries
        )

        # ASSUME: future capacity is distributed to subnational regions based on a
        # linear scaling from the current distribution of capacity
        nuclear_regional_proportion = (
            capacity_current_per_region
            .groupby("country_code")
            .transform(lambda x: x / x.sum())
        )

        capacity_per_region = (
            future_capacity
            .mul(nuclear_regional_proportion, level='country_code', axis=0)
        )

        # Check that we haven't lost any capacity on regionalisation
        assert capacity_per_region.sum(level="country_code").equals(future_capacity)

    capacity_per_region.dropna(how="all").sum(level="id").to_csv(path_to_output)


def _iso3(country_name):
    return pycountry.countries.lookup(country_name).alpha_3


def _get_future_capacity_from_config_file(nuclear_capacity_scenario, resolution, countries):

    nuclear_scenario_df = (
        pd.read_csv(nuclear_capacity_scenario, index_col="country")
        .loc[:, ["min", "max"]]
        .rename(columns=lambda x: f"installed_capacity_nuclear_{x}_MW")
        .reindex(countries)
        .dropna()
    )

    if resolution == "continental":
        # sum over only those countries which are relevant to this workflow run
        future_capacity = nuclear_scenario_df.groupby(lambda x: "EUR").sum()
    else:
        future_capacity = nuclear_scenario_df.rename(index=_iso3)

    return future_capacity.rename_axis(index="country_code").astype(float)


if __name__ == "__main__":
    regionalise_nuclear_capacity(
        path_to_power_plant_database=snakemake.input.power_plant_database,
        path_to_units=snakemake.input.units,
        nuclear_capacity_scenario=snakemake.params.nuclear_capacity_scenario,
        countries=snakemake.params.countries,
        resolution=snakemake.wildcards.resolution,
        path_to_output=snakemake.output[0]
    )
