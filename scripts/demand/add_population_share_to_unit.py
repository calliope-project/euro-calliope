import geopandas as gpd

DRIVER = "GeoJSON"


def add_share(path_to_units_with_population: str, path_to_output: str) -> None:

    """
    Computes a unit's population's share of the national population and adds share to unit.
    For national and continental resolution the population share is 1.0 for each unit.

    Parameters:
        path_to_units_with_population (str):
            Location of geojson file with units
        path_to_output (str):
            Location to store generated geojson file with population shares
    """

    units = gpd.read_file(path_to_units_with_population).set_index("id")
    units["population_share"] = units.groupby("country_code").population_sum.transform(lambda x: x / x.sum())
    units.to_file(path_to_output, driver=DRIVER)


if __name__ == "__main__":
    add_share(
        path_to_units_with_population=snakemake.input.units,
        path_to_output=snakemake.output[0]
    )
