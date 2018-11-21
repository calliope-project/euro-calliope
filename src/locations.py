"""Creates Calliope location files."""
from pathlib import Path
from collections import namedtuple

import pandas as pd
from calliope import AttrDict



SubLocation = namedtuple("SubLocation", ["name", "techs", "eligibilities"])
"""A sub location of an administrative unit.

SubLocations are used to define sizes of areas eligible for renewable electricity
installations.
"""

SUB_LOCATIONS = [
    SubLocation(
        "wind_onshore",
        {"wind_onshore": {"constraints": {"resource": "file=capacityfactors-wind-onshore.csv:{column}"}}},
        ["eligibility_onshore_wind_km2"]
    ),
    SubLocation(
        "wind_offshore",
        {"wind_offshore": {"constraints": {"resource": "file=capacityfactors-wind-offshore.csv:{column}"}}},
        ["eligibility_offshore_wind_km2"]
    ),
    SubLocation(
        "pv_or_wind_farm",
        {"open_field_pv": {"constraints": {"resource": "file=capacityfactors-open-field-pv.csv:{column}"}},
         "wind_onshore": {"constraints": {"resource": "file=capacityfactors-wind-onshore.csv:{column}"}}},
        ["eligibility_onshore_wind_and_pv_km2"]
    ),
    SubLocation(
        "roof_mounted_pv",
        {"roof_mounted_pv": {"constraints": {"resource": "file=capacityfactors-rooftop-pv.csv:{column}"}}},
        ["eligibility_rooftop_pv_km2"]
    )
]


def generate_locations(path_to_land_eligibility_km2, path_to_result):
    """Generate a file that represents locations in Calliope."""
    land_eligibility_km2 = pd.read_csv(path_to_land_eligibility_km2, index_col=0)
    ids = land_eligibility_km2.index

    locations = generate_main_locations(ids)
    locations.union(generate_sub_locations(ids, land_eligibility_km2))
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


def generate_sub_locations(ids, land_eligibility_km2):
    config = AttrDict()
    for id in ids:
        for sub_location in SUB_LOCATIONS:
            available_area = land_eligibility_km2.loc[id, sub_location.eligibilities].sum()
            sub_id = id + "_" + sub_location.name
            sub_location_config = AttrDict.from_yaml_string(
                f"""
                locations:
                    {sub_id}:
                        available_area: {available_area}
                        techs:
                """
            )
            sub_location_config.set_key(
                f"locations.{sub_id}.techs",
                {tech: _replace_column_name(config, id)
                 for tech, config in sub_location.techs.items()}
            )
            config.union(sub_location_config)
    return config


def _replace_column_name(config, super_id):
    config["constraints"]["resource"] = config["constraints"]["resource"].format(column=super_id)
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
    generate_locations(snakemake.input.land_eligibility_km2, snakemake.output[0])
