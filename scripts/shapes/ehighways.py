"""Remixes NUTS, LAU, GADM, and E-Highways shapes to form the units of the analysis."""

import geopandas as gpd
import pandas as pd
import shapely

OUTPUT_DRIVER = "GPKG"


def create_ehighways_shapes(
    path_to_nuts: str,
    path_to_unit_mapping: str,
    nuts_year: int,
    path_to_output: str,
):
    """Aggregate NUTS shapes into ehighways shapes.

    Args:
        path_to_nuts (str): CSV file that specifies NUTS shape data.
        path_to_unit_mapping (str): CSV file containing a mapping between NUTS shapes and ehighways shapes.
        nuts_year (int): The NUTS reference year used in `path_to_nuts_units`.
        path_to_output (str): Path to which ehighways shapes will be saved.

    Assumes:
    - the NUTS data being mapped is NUTS3
    """

    unit_mapping_df = pd.read_csv(path_to_unit_mapping, header=0)

    nuts_year_string = f"NUTS3_{nuts_year}"

    # locations is a series indexed by the NUTS unit with values of the id (e.g., the ehighways unit)
    unit_mapping_df_no_nan = (
        unit_mapping_df.dropna(subset=[nuts_year_string]).set_index(nuts_year_string).id
    )

    nuts_units = gpd.read_file(path_to_nuts, layer="nuts3").set_index("id")
    nuts_units["ehighways_id"] = unit_mapping_df_no_nan.reindex(nuts_units.index)

    ehighways_units = nuts_units.dissolve("ehighways_id")
    ehighways_units.geometry = ehighways_units.geometry.map(_to_multi_polygon)
    ehighways_units = ehighways_units.rename_axis(index="id").assign(
        name=ehighways_units.index, type="ehighways"
    )
    ehighways_units.to_file(path_to_output, driver=OUTPUT_DRIVER, layer="ehighways")


def _to_multi_polygon(geometry):
    "Handles result of dissolving geometries that don't share a border by converting them to multipolygons"
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.polygon.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


if __name__ == "__main__":
    create_ehighways_shapes(
        path_to_nuts=snakemake.input.nuts,
        path_to_unit_mapping=snakemake.input.nuts_to_ehighways_units,
        nuts_year=snakemake.params.nuts_year,
        path_to_output=snakemake.output[0],
    )
