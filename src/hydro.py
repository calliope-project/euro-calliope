import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import pycountry

WGS_84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


def main(path_to_plants, path_to_locations, path_to_phs_storage_capacities, path_to_output):
    locations = gpd.read_file(path_to_locations).to_crs(WGS_84).set_index("id")
    plants = pd.read_csv(path_to_plants, index_col="id")

    hror = capacities_per_location(plants[plants.type == "HROR"].copy(), locations, tech_type="hror",
                                   fill_storage_capacity=False)
    hdam = capacities_per_location(plants[plants.type == "HDAM"].copy(), locations, tech_type="hdam")
    hphs = capacities_per_location(plants[plants.type == "HPHS"].copy(), locations, tech_type="hphs")
    hphs["storage_capacity_hphs_MWh"] = scale_phs_storage_capacities(
        hphs=hphs,
        locations=locations,
        national_storage_capacities=read_national_phs_storage_capacities(path_to_phs_storage_capacities, locations)
    )

    pd.concat([hror, hdam, hphs], axis="columns").to_csv(
        path_to_output,
        header=True,
        index=True
    )


def capacities_per_location(plants, locations, tech_type, fill_storage_capacity=True):
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


def read_national_phs_storage_capacities(path_to_data, locations):
    data = pd.read_csv(path_to_data, index_col=0)
    data.index = [pycountry.countries.lookup(iso2).alpha_3 for iso2 in data.index]
    data["storage-capacity-mwh"] = data["storage-capacity-gwh"] * 1000
    if (len(locations.index) == 1) and (locations.index[0] == "EUR"): # special case for continental level
        data = pd.DataFrame(index=["EUR"], data=data.sum(axis=0).to_dict())
    return data.reindex(locations.country_code.unique(), fill_value=0) # FIXME kills capacity in Romania


def scale_phs_storage_capacities(hphs, locations, national_storage_capacities):
    """Scale PHS storage capacities to match (Geth et al., 2015).

    Storage capacities of pumped hydro within the JRC data base seem too high. Thus, we are scaling
    them here, so that the national numbers of (Geth et al., 2015) are fulfilled.

    Geth, F., Brijs, T., Kathan, J., Driesen, J., Belmans, R., 2015. An overview of large-scale
    stationary electricity storage plants in Europe: Current status and new developments. Renewable
    and Sustainable Energy Reviews 52, 1212–1227. https://doi.org/10.1016/j.rser.2015.07.145
    """
    locations = pd.concat([
        locations,
        hphs.groupby(locations.country_code)
            .storage_capacity_hphs_MWh
            .transform(lambda x: x / x.sum())
            .rename("storage_capacity_share")
    ], axis=1)
    locations = locations.merge(national_storage_capacities, right_index=True, left_on="country_code")
    new_storage_capacities = locations["storage_capacity_share"] * locations["storage-capacity-mwh"]
    return new_storage_capacities.rename("storage_capacity_hphs_MWh").fillna(0)


if __name__ == "__main__":
    main(
        path_to_plants=snakemake.input.plants,
        path_to_locations=snakemake.input.locations,
        path_to_phs_storage_capacities=snakemake.input.phs_storage_capacities,
        path_to_output=snakemake.output[0]
    )
