"""Creates Calliope location files."""
from pathlib import Path
from collections import namedtuple

import pandas as pd
import geopandas as gpd
from calliope import AttrDict


SubLocation = namedtuple("SubLocation", ["name", "techs", "eligibilities"])
"""A sub location of an administrative unit.

SubLocations are used to define sizes of areas eligible for renewable electricity
installations.
"""

SUB_LOCATIONS = [
    SubLocation(
        "wind_onshore",
        ["wind_onshore"],
        ["eligibility_onshore_wind_other_km2", "eligibility_onshore_wind_farmland_km2",
         "eligibility_onshore_wind_forest_km2", "eligibility_onshore_wind_other_protected_km2",
         "eligibility_onshore_wind_farmland_protected_km2", "eligibility_onshore_wind_forest_protected_km2"]
    ),
    SubLocation(
        "pv_or_wind_farm",
        ["pv", "wind_onshore"],
        ["eligibility_onshore_wind_and_pv_other_km2", "eligibility_onshore_wind_and_pv_farmland_km2",
         "eligibility_onshore_wind_and_pv_other_protected_km2",
         "eligibility_onshore_wind_and_pv_farmland_protected_km2"]
    ),
    SubLocation(
        "roof_mounted_pv",
        ["pv"],
        ["eligibility_rooftop_pv_km2"]
    )
]


def generate_locations(path_to_ids, path_to_eligibility, path_to_result):
    """Generate a file that represents locations in Calliope.

    Takes as input a GeoJSON that defines administrative units and creates one
    Calliope location for each administrative unit.
    """
    units = gpd.read_file(path_to_ids)
    land_eligibility = pd.read_csv(path_to_eligibility, index_col=0)
    ids = units.id.values

    locations = generate_main_locations(ids)
    locations.union(generate_sub_locations(ids))
    locations.union(generate_sub_location_area_constraints(ids, land_eligibility))
    locations.union(generate_sub_links(ids))
    locations.to_yaml(Path(path_to_result))


def generate_main_locations(ids):
    ids = ','.join(ids)
    return AttrDict.from_yaml_string(
        f"""
        locations:
            {ids}:
                techs:
                    demand_elec:
                    battery:
        """
    )


def generate_sub_locations(ids):
    config = AttrDict()
    for sub_location in SUB_LOCATIONS:
        sub_ids = [id + "_" + sub_location.name for id in ids]
        sub_ids = ','.join(sub_ids)
        sub_location_config = AttrDict.from_yaml_string(
            f"""
            locations:
                {sub_ids}:
                    techs:
            """
        )
        sub_location_config.set_key(f"locations.{sub_ids}.techs", {tech: None for tech in sub_location.techs})
        config.union(sub_location_config)
    return config


def generate_sub_location_area_constraints(ids, land_eligibility):
    config = AttrDict()
    for id in ids:
        for sub_location in SUB_LOCATIONS:
            available_area = land_eligibility.loc[id, sub_location.eligibilities].sum()
            sub_id = id + "_" + sub_location.name
            config.union(
                AttrDict.from_yaml_string(
                    f"""
                    locations:
                        {sub_id}.available_area: {available_area}
                    """
                )
            )
    return config


def generate_sub_links(ids):
    config = AttrDict()
    for id in ids:
        for sub_location in SUB_LOCATIONS:
            sub_id = id + "_" + sub_location.name
            config.union(
                AttrDict.from_yaml_string(
                    f"""
                    links:
                        {id},{sub_id}.techs:
                            free_transmission:
                    """
                )
            )
    return config


if __name__ == "__main__":
    generate_locations(snakemake.input.ids, snakemake.input.land_eligibility, snakemake.output[0])
