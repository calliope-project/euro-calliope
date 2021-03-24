import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

WGS_84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


def main(path_to_plants, path_to_locations, path_to_output):
    locations = gpd.read_file(path_to_locations).to_crs(WGS_84).set_index("id")
    plants = pd.read_csv(path_to_plants, index_col="id")

    hror = capacities_per_location(plants[plants.type == "HROR"].copy(), locations, tech_type="hror",
                                   fill_storage_capacity=False)
    hdam = capacities_per_location(plants[plants.type == "HDAM"].copy(), locations, tech_type="hdam")
    hphs = capacities_per_location(plants[plants.type == "HPHS"].copy(), locations, tech_type="hphs")

    pd.concat([hror, hdam, hphs], axis="columns").to_csv(
        path_to_output,
        header=True,
        index=True
    )


def capacities_per_location(plants, locations, tech_type, fill_storage_capacity=True):
    if fill_storage_capacity:
        plants = fill_missing_storage_capacity_values(plants.copy())
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
                  .rename(columns={"installed_capacity_MW": f"installed_capacity_{tech_type}_MW"})
                  .rename(columns={"storage_capacity_MWh": f"storage_capacity_{tech_type}_MWh"}))


def fill_missing_storage_capacity_values(plants):
    # ASSUME country specific median E/P ratio for missing values, global median where no country specific available
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
        path_to_output=snakemake.output[0]
    )
