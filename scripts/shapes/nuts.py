"""Preprocessing of raw NUTS data to bring it into normalised form."""

import logging
from typing import Union

import fiona
import geopandas as gpd
import numpy as np
import shapely.geometry
from eurocalliopelib import utils

OUTPUT_DRIVER = "GPKG"

LOGGER = logging.getLogger("nuts.py")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(name)s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)


def clean_nuts(
    path_to_nuts: str,
    path_to_output: str,
    crs: str,
    study_area: shapely.geometry.Polygon,
    all_countries: list[str],
    schema: dict,
):
    """Preprocess downloaded NUTS data into a format that aligns with other shapefile data sources (e.g. GADM).

    Args:
        path_to_nuts (str): Path to downloaded NUTS data.
        path_to_output (str): Path to which output GeoPackage will be saved.
        crs (str): Study CRS. Output data will be stored in this CRS.
        study_area (shapely.geometry.Polygon): Study area to which NUTS data will be clipped.
        all_countries (list[str]): List of countries to which NUTS data will be clipped.
        schema (dict): Schema to which the processed data should align.
    """
    gdf_nuts = gpd.read_file(path_to_nuts)

    gdf_nuts_clean = (
        gpd.GeoDataFrame(
            geometry=gdf_nuts.geometry.apply(_to_multi_polygon),
            data={
                "id": gdf_nuts.NUTS_ID,
                "name": gdf_nuts.NAME_LATN,
                "type": np.where(gdf_nuts.LEVL_CODE == 0, "country", np.nan),
                "layer": "nuts" + gdf_nuts.LEVL_CODE.astype(str),
                "proper": True,
                "country_code": gdf_nuts.CNTR_CODE.apply(utils.eu_country_code_to_iso3),
            },
            crs=gdf_nuts.crs,
        )
        .to_crs(crs)
        .pipe(_fix_country_names)
        .pipe(_clip_nuts, study_area, all_countries)
    )

    gdf_nuts_clean.groupby("layer").apply(
        _to_file, path_to_output=path_to_output, schema=schema
    )

    _test_id_uniqueness(path_to_output)


def _clip_nuts(
    gdf_nuts: gpd.GeoDataFrame,
    study_area: shapely.geometry.Polygon,
    all_countries: list[str],
) -> gpd.GeoDataFrame:
    """Keep only geometries that are in the study area.

    Remove:
    * countries not in study country codes
    * entire shapes is not in study area bounding box
    * parts of multipolygons that are not in the study area (e.g. remove Caribbean islands from the France country shape)

    """
    countries = utils.convert_valid_countries(all_countries, output="alpha3").values()
    to_keep = gdf_nuts.intersects(study_area) & gdf_nuts.country_code.isin(countries)
    to_drop = gdf_nuts[~to_keep]
    LOGGER.info(f"Removing NUTS regions outside study area: {to_drop.id.values}")
    gdf_nuts_clipped = gdf_nuts[to_keep].apply(
        _all_parts_in_study_area, study_area=study_area, axis=1
    )
    return gdf_nuts_clipped


def _fix_country_names(gdf_nuts: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    "Use pycountry names and ISO3 ID for NUTS0 shapes, instead of those supplied by NUTS"

    gdf_nuts["name"] = np.where(
        gdf_nuts.layer == "nuts0",
        gdf_nuts.country_code.apply(utils.convert_country_code, output="name"),
        gdf_nuts["name"],
    )
    gdf_nuts["id"] = np.where(
        gdf_nuts.layer == "nuts0", gdf_nuts.country_code, gdf_nuts["id"]
    )
    return gdf_nuts


def _all_parts_in_study_area(
    feature: gpd.GeoSeries, study_area: shapely.geometry.Polygon
) -> gpd.GeoSeries:
    unit = feature.geometry
    if not study_area.contains(unit):
        LOGGER.info(f"Removing parts of {feature.id} outside of study area.")
        new_unit = shapely.geometry.MultiPolygon([
            polygon for polygon in unit.geoms if study_area.contains(polygon)
        ])
        feature["geometry"] = new_unit
    return feature


def _to_multi_polygon(
    geometry: Union[dict, shapely.geometry.Polygon],
) -> shapely.geometry.MultiPolygon:
    "Convert all input shapes into MultiPolygons"
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


def _to_file(gdf: gpd.GeoDataFrame, path_to_output: str, schema: dict):
    "Save a layer to a geopackage"
    gdf.drop("layer", axis=1).to_file(
        path_to_output, driver=OUTPUT_DRIVER, layer=gdf.name, schema=schema
    )


def _test_id_uniqueness(path_to_file: str):
    for layer_name in fiona.listlayers(path_to_file):
        assert not gpd.read_file(path_to_file, layer=layer_name).id.duplicated().any()


if __name__ == "__main__":
    clean_nuts(
        path_to_nuts=snakemake.input.zipped,
        path_to_output=snakemake.output[0],
        crs=snakemake.params.crs,
        schema=snakemake.params.schema,
        all_countries=snakemake.params.all_countries,
        study_area=shapely.geometry.box(
            minx=snakemake.params.x_min,
            maxx=snakemake.params.x_max,
            miny=snakemake.params.y_min,
            maxy=snakemake.params.y_max,
        ),
    )
