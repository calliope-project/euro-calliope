"""Remixes NUTS, LAU, and GADM data to form the units of the analysis."""

import fiona
import geopandas as gpd
import pandas as pd
import pycountry
import shapely
from eurocalliopelib import utils

DRIVER = "GeoJSON"


def remix_units(
    path_to_nuts: str,
    path_to_gadm: str,
    path_to_output: str,
    resolution,
    layer_config,
    all_countries,
    path_to_statistical_to_custom_units,
):
    """Remixes NUTS, LAU, and GADM data to form the units of the analysis."""
    # TODO: do I need to understand the validation functions?
    source_layers = _read_source_layers(path_to_nuts, path_to_gadm)
    _validate_source_layers(source_layers)
    _validate_layer_config(all_countries, layer_config, resolution)
    units = _build_layer(layer_config, source_layers)
    _validate_layer(units, resolution, all_countries)
    if resolution == "continental":
        # sum all of the units together into one block
        units = _continental_layer(units)
    elif path_to_statistical_to_custom_units != []:
        # this code is only executed if there is a [statistical-to-custom-units file corresponding to the chosen resolution] in the config
        # sum units together according to custom mapping
        units = _custom_layer(units, path_to_statistical_to_custom_units)
    _write_layer(units, path_to_output)


def _custom_layer(base_units, path_to_statistical_to_custom_units):
    """Combine the lower-level statistical units (aspecified in layer_config)
    into the custom units, according to statistical_to_custom_units, using _merge_units
    function.

    Inputs:
    - path_to_statistical_to_custom_units: csv file that specifies which NUTS and GADM
        units contribute to each custom unit
    - layer: the GeoDataFrame covering the whole scope, where each country is formed
        from units corresponding to the resolution specified in layer_config
    Outputs:
    -

    Assumes:
    - the base-level statistical units in layer and are either NUTS3 or whole countries (NUTS0 or GADM0)
    - any units which are whole countries will also be whole countries at the custom resolution
    """

    stat_to_custom_units_df = pd.read_csv(path_to_statistical_to_custom_units, header=0)
    nuts_year = 2013  # This was a variable in Bryn's SC EC, but hardcoded to 2013 in EC. Consider that the year needs to match for both the (downloaded nuts units) and (nuts to custom region mapping)
    nuts_year_string = f"NUTS3_{nuts_year}"

    # locations is a series indexed by the statistical unit (NUTS or GADM) with values of the id (e.g., the ehighways unit)
    locations = (
        stat_to_custom_units_df.dropna(subset=[nuts_year_string])
        .set_index(
            stat_to_custom_units_df.dropna(subset=[nuts_year_string])[nuts_year_string]
        )
        .id.append(
            stat_to_custom_units_df[stat_to_custom_units_df.Source != "NUTS3"]
            .set_index("country_code")
            .rename(index=utils.eu_country_code_to_iso3)
            .id
        )
    )

    units = _merge_units(base_units, locations)
    return units


def _merge_units(base_units, locations):
    """Summary:
    - uses 'dissolve' method of GeoPandas to combine base units into custom units
    - _to_multi_polygon function handles shapes with same unit but no shared border
    """

    def _to_multi_polygon(geometry):
        if isinstance(geometry, dict):
            geometry = shapely.geometry.shape(geometry)
        if isinstance(geometry, shapely.geometry.polygon.Polygon):
            return shapely.geometry.MultiPolygon(polygons=[geometry])
        else:
            return geometry

    base_units = base_units.set_index("id")
    base_units["custom_unit"] = (locations.reindex(base_units.index)).dropna()
    base_units.loc[~base_units.is_valid, "geometry"] = base_units.loc[
        ~base_units.is_valid, "geometry"
    ].buffer(0)
    units = base_units.dissolve("custom_unit")
    units.geometry = units.geometry.map(_to_multi_polygon)
    units.index.rename("id", inplace=True)
    units = units.reset_index()
    units.loc[units["type"].isnull(), "name"] = "custom_unit"
    units["type"].fillna("custom_unit", inplace=True)

    return units


def _build_layer(layer_config, source_layers):
    """
    Inputs:
    - layer_config: a mapping of which statistical unit resolution each country should use as its base
    - source_layers: a dict with keys as each statistical unit type (nuts[0-3] or gadmn[0-3])
        and values as a GeoDataFrame representing the collection of those units
    Returns:
    - layer: the GeoDataFrame covering the whole scope in source_layers, where each
        country is formed from units corresponding to the resolution specified in
        layer_config
    """
    crs = [layer.crs for layer in source_layers.values()][0]
    layer = pd.concat([
        source_layers[source_layer][
            source_layers[source_layer].country_code == _iso3(country)
        ]
        for country, source_layer in layer_config.items()
    ])
    assert isinstance(layer, pd.DataFrame)
    return gpd.GeoDataFrame(layer, crs=crs)


def _read_source_layers(path_to_nuts, path_to_gadm):
    """Returns:
    - source_layers: a dict with keys as each statistical unit type (nuts[0-3] or gadmn[0-3])
        and values as a geodataframe representing the collection of those units.
    """
    source_layers = {
        layer_name: gpd.read_file(path_to_nuts, layer=layer_name)
        for layer_name in fiona.listlayers(path_to_nuts)
    }
    source_layers.update({
        layer_name: gpd.read_file(path_to_gadm, layer=layer_name)
        for layer_name in fiona.listlayers(path_to_gadm)
    })
    return source_layers


def _validate_source_layers(source_layers):
    crs = [layer.crs for layer in source_layers.values()]
    assert not crs or crs.count(crs[0]) == len(
        crs
    ), "Source layers have different crs. They must match."


def _validate_layer_config(all_countries, layer_config, layer_name):
    assert all(country in layer_config for country in all_countries), (
        f"Layer {layer_name} is not correctly " "defined."
    )


def _validate_layer(layer, layer_name, countries):
    assert all(
        _iso3(country) in layer.country_code.unique() for country in countries
    ), f"Countries are missing in layer {layer_name}."


def _iso3(country_name):
    return pycountry.countries.lookup(country_name).alpha_3


def _continental_layer(layer):
    # special case all Europe
    layer = layer.dissolve(by=[1 for idx in layer.index])
    index = layer.index[0]
    layer.loc[index, "id"] = "EUR"
    layer.loc[index, "country_code"] = "EUR"
    layer.loc[index, "name"] = "Europe"
    layer.loc[index, "type"] = "continent"
    layer.loc[index, "proper"] = 1
    return layer


def _write_layer(gdf, path_to_file):
    gdf.to_file(path_to_file, driver=DRIVER)


if __name__ == "__main__":
    resolution = snakemake.wildcards[0]
    remix_units(
        path_to_nuts=snakemake.input.nuts,
        path_to_gadm=snakemake.input.gadm,
        path_to_output=snakemake.output[0],
        all_countries=snakemake.params.all_countries,
        layer_config=snakemake.params.layer_configs[resolution],
        resolution=resolution,
        path_to_statistical_to_custom_units=snakemake.input.statistical_to_custom_units,
    )
