import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import jinja2

WGS_84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

TEMPLATE = """
locations:
    {% for location_index in hror.index %}
    {{ location_index }}:
        techs:
            hydro_run_of_river:
                constraints:
                    energy_cap_equals: {{ hror.loc[location_index, "installed_capacity_MW"] * scaling_factor }} # [{{ 1 / scaling_factor }} MW]
            hydro_reservoir:
                constraints:
                    energy_cap_equals: {{ hdam.loc[location_index, "installed_capacity_MW"] * scaling_factor }} # [{{ 1 / scaling_factor }} MW]
                    storage_cap_equals: {{ hdam.loc[location_index, "storage_capacity_MWh"] * scaling_factor }} # [{{ 1 / scaling_factor }} MWh]

            pumped_hydro:
                constraints:
                    energy_cap_equals: {{ hphs.loc[location_index, "installed_capacity_MW"] * scaling_factor }} # [{{ 1 / scaling_factor }} MW]
                    storage_cap_equals: {{ hphs.loc[location_index, "storage_capacity_MWh"] * scaling_factor }} # [{{ 1 / scaling_factor }} MWh]
    {% endfor %}
"""


def main(path_to_plants, path_to_locations, scaling_factor, path_to_output):
    locations = gpd.read_file(path_to_locations).to_crs(WGS_84).set_index("id")
    plants = pd.read_csv(path_to_plants, index_col="id")

    yaml_capacities = jinja2.Template(TEMPLATE).render(
        hror=capacities_per_location(plants[plants.type == "HROR"].copy(), locations, fill_storage_capacity=False),
        hdam=capacities_per_location(plants[plants.type == "HDAM"].copy(), locations),
        hphs=capacities_per_location(plants[plants.type == "HPHS"].copy(), locations),
        scaling_factor=scaling_factor
    )
    with open(path_to_output, "w") as result_file:
        result_file.write(yaml_capacities)


def capacities_per_location(plants, locations, fill_storage_capacity=True):
    plants = fill_missing_storage_capacity_values(plants)
    plant_centroids = gpd.GeoDataFrame(
        crs=WGS_84,
        geometry=list(map(Point, zip(plants.lon, plants.lat))),
        index=plants.index
    )
    location_of_plant = gpd.sjoin(plant_centroids, locations, how="left", op='intersects')["index_right"]
    location_of_plant = location_of_plant[~location_of_plant.index.duplicated()]
    return (plants.groupby(location_of_plant)
                  .agg({"installed_capacity_MW": sum, "storage_capacity_MWh": sum})
                  .reindex(index=locations.index, fill_value=0)
                  .rename(index=lambda idx: idx.replace(".", "-")))


def fill_missing_storage_capacity_values(plants):
    # ASSUME country specific median E/P ratio for missing values, global median where no country specific available
    # TODO revisit this ^ assumption
    e_to_p = plants.storage_capacity_MWh / plants.installed_capacity_MW
    based_on_global = e_to_p.median() * plants.installed_capacity_MW
    based_on_country_specific = plants.merge(
        e_to_p.groupby(plants.country_code).median().rename("country_specific"),
        left_on="country_code",
        right_index=True,
        how="left"
    ).loc[:, "country_specific"] * plants.installed_capacity_MW
    plants["storage_capacity_MWh"] = (plants["storage_capacity_MWh"].where(pd.notnull, other=based_on_country_specific)
                                                                    .where(pd.notnull, other=based_on_global))
    assert not plants["storage_capacity_MWh"].isnull().any()
    return plants


if __name__ == "__main__":
    main(
        path_to_plants=snakemake.input.plants,
        path_to_locations=snakemake.input.locations,
        scaling_factor=snakemake.params.scaling_factor,
        path_to_output=snakemake.output[0]
    )
