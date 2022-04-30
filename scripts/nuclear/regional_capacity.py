import pandas as pd
import shapely.geometry
import geopandas as gpd
import pycountry


def regionalise_nuclear_capacity(
    path_to_power_plant_database, path_to_units,
    nuclear_scenario_config, path_to_output
):
    """
    Use current geolocations of nuclear capacity in Europe to assign nuclear
    capacity to Euro-Calliope regions.

    Args:
        path_to_power_plant_database (str):
            location of JRC powerplant database "OPEN_UNITS" CSV file
        path_to_units (str):
            location of euro-calliope regions geojson/shapefile
        nuclear_scenario_config (dict):
            nuclear scenario name and (if scenario != 'current')
            minimum and maximum capacities of future nuclear capacity per country in MW
            (e.g. {'scenario: 'future_capacity', 'national_capacities': {'future_capacity': {'France': {'min': 1000, 'max': 5000}}}})
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
    scenario = nuclear_scenario_config["scenario"]
    if scenario == "current":
        capacity_per_region = (
            capacity_current_per_region
            .to_frame("installed_capacity_nuclear_equals_MW")
        )
    else:
        nuclear_regional_proportion = capacity_current_per_region.div(
            capacity_current_per_region.sum(level='country_code')
        )
        future_capacity = (
            pd.DataFrame(nuclear_scenario_config["national_capacities"][scenario])
            .transpose()
            .rename(index=_iso3, columns=lambda x: f"installed_capacity_nuclear_{x}_MW")
            .rename_axis(index="country_code")
        )

        capacity_per_region = (
            future_capacity
            .mul(nuclear_regional_proportion, level='country_code', axis=0)
        )

    capacity_per_region.dropna(how="all").sum(level="id").to_csv(path_to_output)


def _iso3(country_name):
    return pycountry.countries.lookup(country_name).alpha_3


if __name__ == "__main__":
    regionalise_nuclear_capacity(
        path_to_power_plant_database=snakemake.input.power_plant_database,
        path_to_units=snakemake.input.units,
        nuclear_scenario_config=snakemake.params.nuclear_scenario_config,
        path_to_output=snakemake.output[0]
    )
