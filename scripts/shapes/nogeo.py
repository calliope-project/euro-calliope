import pandas as pd
import geopandas as gpd


def remove_geo_information(path_to_units, path_to_units_without_geo_information):
    units = gpd.read_file(path_to_units).set_index("id")
    units = pd.DataFrame(units)
    del units["geometry"]
    units.to_csv(path_to_units_without_geo_information, index=True, header=True)


if __name__ == "__main__":
    remove_geo_information(
        path_to_units=snakemake.input.units,
        path_to_units_without_geo_information=snakemake.output[0]
    )
