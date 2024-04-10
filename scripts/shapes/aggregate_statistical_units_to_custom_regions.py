import geopandas as gpd
import pandas as pd
import shapely.geometry
from eurocalliopelib import utils

DRIVER = "GeoJSON"


def dissolve_nuts3(
    path_to_disaggregated_units, path_to_output, nuts_year, path_to_nuts_to_regions
):
    """Dissolve NUTS3 data to custom regions."""

    nuts_to_regions_df = pd.read_csv(path_to_nuts_to_regions, header=0)
    # Any clusters which are whole countries will be whole countries at the subregional scale too
    nuts_year_string = f"NUTS3_{nuts_year}"
    locations = (
        nuts_to_regions_df.dropna(subset=[nuts_year_string])
        .set_index(
            nuts_to_regions_df.dropna(subset=[nuts_year_string])[nuts_year_string]
        )
        .id.append(
            nuts_to_regions_df[nuts_to_regions_df.Source != "NUTS3"]
            .set_index("country_code")
            .rename(index=utils.eu_country_code_to_iso3)
            .id
        )
    )

    shapes = _update_units(path_to_disaggregated_units, locations)
    shapes.to_file(path_to_output, driver=DRIVER)


def _update_units(path_to_units, locations):
    def _to_multi_polygon(geometry):
        if isinstance(geometry, dict):
            geometry = shapely.geometry.shape(geometry)
        if isinstance(geometry, shapely.geometry.polygon.Polygon):
            return shapely.geometry.MultiPolygon(polygons=[geometry])
        else:
            return geometry

    units = gpd.read_file(path_to_units).set_index("id")
    units["custom_region"] = (locations.reindex(units.index)).dropna()
    units.loc[~units.is_valid, "geometry"] = units.loc[
        ~units.is_valid, "geometry"
    ].buffer(0)
    units = units.dissolve("custom_region")
    units.geometry = units.geometry.map(_to_multi_polygon)
    units.index.rename("id", inplace=True)
    units = units.reset_index()
    units.loc[units["type"].isnull(), "name"] = "custom_region"
    units["type"].fillna("custom_region", inplace=True)

    return units


if __name__ == "__main__":
    dissolve_nuts3(
        path_to_disaggregated_units=snakemake.input.disaggregated_units,
        path_to_nuts_to_regions=snakemake.input.nuts_to_regions,
        nuts_year=int(snakemake.params.nuts_year),
        path_to_output=snakemake.output[0],
    )
