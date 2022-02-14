"""Creates Calliope location files."""
import geopandas as gpd

from eurocalliopelib.template import parametrise_template
from eurocalliopelib.geo import WGS84, EPSG3035


def construct_locations(
    path_to_shapes, path_to_template, resolution, path_to_output_yaml, path_to_output_csv
):
    """Generate a file that represents locations in Calliope."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id")
    )
    locations = (
        locations
        .assign(centroid=locations.to_crs(EPSG3035).centroid.to_crs(WGS84).rename("centroid"))
        .loc[:, ["name", "centroid"]]
    )

    locations = locations.assign(id=locations.index.str.replace(".", "-", regex=False)).set_index("id")

    parametrise_template(
        path_to_template, path_to_output_yaml,
        locations=locations, resolution=resolution
    )

    locations.name.to_csv(path_to_output_csv, index=True, header=True)


if __name__ == "__main__":
    construct_locations(
        path_to_shapes=snakemake.input.shapes,
        path_to_template=snakemake.input.template,
        resolution=snakemake.wildcards.resolution,
        path_to_output_yaml=snakemake.output.yaml,
        path_to_output_csv=snakemake.output.csv,
    )
