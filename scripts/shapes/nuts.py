"""Preprocessing of raw NUTS data to bring it into normalised form."""
import zipfile

import fiona
import fiona.transform
import shapely.geometry
import geopandas as gpd
import pandas as pd
import pycountry

from eurocalliopelib import utils


OUTPUT_DRIVER = "GPKG"
LAYER_NAME = "nuts{layer_id}"


def merge(path_to_shapes, path_to_attributes, path_to_output):
    """Merge NUTS shapes with attributes."""
    shapes = gpd.read_file(path_to_shapes)
    shapes.geometry = shapes.geometry.map(_to_multi_polygon)
    attributes = gpd.read_file(path_to_attributes)
    attributes = pd.DataFrame(attributes) # to be able to remove geo information
    del attributes["geometry"]
    shapes.merge(attributes, on="NUTS_ID", how="left").to_file(path_to_output, driver=OUTPUT_DRIVER)


def normalise(path_to_nuts, path_to_output, crs, study_area, all_countries, schema):
    """Normalises raw NUTS data.

    Raw data contains all NUTS layers in one layer of one shapefile. The output
    of this function corresponds to the form the data is used in this analysis,
    where each geographical layer is stored in one layer of a GeoPackage.
    """
    with fiona.open(path_to_nuts, "r") as nuts_file:
        for layer_id in range(4):
            print("Building layer {}...".format(layer_id))
            _write_layer(nuts_file, crs, study_area, all_countries, path_to_output, layer_id, schema)
    _test_id_uniqueness(path_to_output)


def _write_layer(nuts_file, crs, study_area, all_countries, path_to_output, layer_id, schema):
    with fiona.open(path_to_output,
                    "w",
                    crs=crs,
                    schema=schema,
                    driver=OUTPUT_DRIVER,
                    layer=LAYER_NAME.format(layer_id=layer_id)) as result_file:
        result_file.writerecords(_layer_features(nuts_file, crs, study_area, all_countries, layer_id))


def _layer_features(nuts_file, crs, study_area, all_countries, layer_id):
    for feature in filter(_in_layer_and_in_study_area(layer_id, study_area, all_countries), nuts_file):
        new_feature = {}
        new_feature["properties"] = {}
        new_feature["properties"]["country_code"] = utils.eu_country_code_to_iso3(feature["properties"]["NUTS_ID"][:2])
        new_feature["properties"]["id"] = feature["properties"]["NUTS_ID"]
        new_feature["properties"]["name"] = feature["properties"]["NAME_LATN"]
        new_feature["properties"]["type"] = "country" if layer_id == 0 else None
        new_feature["properties"]["proper"] = True
        new_feature["geometry"] = _all_parts_in_study_area_and_crs(feature, nuts_file.crs, crs, study_area)
        if layer_id == 0:
            new_feature = _fix_country_feature(new_feature)
        yield new_feature


def _fix_country_feature(feature):
    # * IDs should have three letters instead of two
    # * many country names are broken or missing
    feature["properties"]["id"] = utils.eu_country_code_to_iso3(feature["properties"]["id"])
    feature["properties"]["name"] = pycountry.countries.lookup(feature["properties"]["id"]).name
    return feature


def _all_parts_in_study_area_and_crs(feature, src_crs, dst_crs, study_area):
    unit = _to_multi_polygon(feature["geometry"])
    if not study_area.contains(unit):
        print("Removing parts of {} outside of study area.".format(feature["properties"]["NUTS_ID"]))
        new_unit = shapely.geometry.MultiPolygon([polygon for polygon in unit.geoms
                                                  if study_area.contains(polygon)])
        unit = new_unit
    geometry = shapely.geometry.mapping(unit)
    return fiona.transform.transform_geom(
        src_crs=src_crs,
        dst_crs=dst_crs,
        geom=geometry
    )


def _in_layer_and_in_study_area(layer_id, study_area, all_countries):
    def _in_layer_and_in_study_area(feature):
        return _in_layer(layer_id, feature) and _in_study_area(study_area, all_countries, feature)
    return _in_layer_and_in_study_area


def _in_layer(layer_id, feature):
    if feature["properties"]["STAT_LEVL_"] == layer_id:
        return True
    else:
        return False


def _in_study_area(study_area, all_countries, feature):
    countries = [pycountry.countries.lookup(country) for country in all_countries]
    unit = shapely.geometry.shape(feature["geometry"])
    country = pycountry.countries.lookup(utils.eu_country_code_to_iso3(feature["properties"]["NUTS_ID"][:2]))
    if (country in countries) and (study_area.contains(unit) or study_area.intersects(unit)):
        return True
    else:
        print("Removing {} as it is outside of study area.".format(feature["properties"]["NUTS_ID"]))
        return False


def _to_multi_polygon(geometry):
    if isinstance(geometry, dict):
        geometry = shapely.geometry.shape(geometry)
    if isinstance(geometry, shapely.geometry.polygon.Polygon):
        return shapely.geometry.MultiPolygon(polygons=[geometry])
    else:
        return geometry


def _test_id_uniqueness(path_to_file):
    for layer_name in fiona.listlayers(path_to_file):
        assert not gpd.read_file(path_to_file, layer=layer_name).id.duplicated().any()


if __name__ == "__main__":
    with zipfile.ZipFile(snakemake.input.zipped, 'r') as zipped:
        zipped.extractall("./build")
    TMP_FILE = "./build/raw-nuts.gpkg"
    merge(
        path_to_shapes="./build/NUTS_2013_01M_SH/data/NUTS_RG_01M_2013.shp",
        path_to_attributes="./build/NUTS_2013_01M_SH/data/NUTS_AT_2013.dbf",
        path_to_output=TMP_FILE
    )
    normalise(
        path_to_nuts=TMP_FILE,
        path_to_output=snakemake.output[0],
        crs=snakemake.params.crs,
        schema=snakemake.params.schema,
        all_countries=snakemake.params.all_countries,
        study_area=shapely.geometry.box(
            minx=snakemake.params.x_min,
            maxx=snakemake.params.x_max,
            miny=snakemake.params.y_min,
            maxy=snakemake.params.y_max
        )
    )
