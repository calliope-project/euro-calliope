"""Creates Calliope location files."""
from pathlib import Path

import geopandas as gpd
from calliope import AttrDict


def generate_locations(path_to_source, path_to_result):
    """Generate a file that represents locations in Calliope.

    Takes as input a GeoJSON that defines administrative units and creates one
    Calliope location for each administrative unit.
    """
    units = gpd.read_file(path_to_source)
    ids = ','.join(units.id.values)
    locations = AttrDict.from_yaml_string(
        f"""
        locations:
            {ids}:
                techs:
                    demand_elec:
                    pv:
                    wind_onshore:
                    battery:
        """
    )
    locations.to_yaml(Path(path_to_result))


if __name__ == "__main__":
    generate_locations(snakemake.input[0], snakemake.output[0])
