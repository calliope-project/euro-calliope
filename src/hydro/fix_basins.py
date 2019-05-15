import geopandas as gpd

DRIVER = "GPKG"


def fix_basins(path_to_basins, path_to_output):
    """Fix the basins shapes which are invalid.

    Following the advice given here:
    https://github.com/Toblerity/Shapely/issues/344
    """
    basins = gpd.read_file(path_to_basins)
    basins.geometry = basins.geometry.map(_buffer_if_necessary)
    basins.to_file(path_to_output, driver=DRIVER)


def _buffer_if_necessary(shape):
    if not shape.is_valid:
        shape = shape.buffer(0.0)
    assert shape.is_valid
    return shape


if __name__ == "__main__":
    fix_basins(
        path_to_basins=snakemake.input.basins,
        path_to_output=snakemake.output[0]
    )
