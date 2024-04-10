"""Remixes NUTS, LAU, and GADM data to form the units of the analysis."""

import fiona
import geopandas as gpd
import pandas as pd
import pycountry

DRIVER = "GeoJSON"


def remix_units(
    path_to_nuts, path_to_gadm, path_to_output, resolution, layer_config, all_countries
):
    """Remixes NUTS, LAU, and GADM data to form the units of the analysis."""
    source_layers = _read_source_layers(path_to_nuts, path_to_gadm)
    _validate_source_layers(source_layers)
    _validate_layer_config(all_countries, layer_config, resolution)
    layer = _build_layer(layer_config, source_layers)
    _validate_layer(layer, resolution, all_countries)
    if resolution == "continental":  # treat special case
        layer = _continental_layer(layer)
    elif mapping is not None:
        layer = _custom_layer(layer, mapping)
    _write_layer(layer, path_to_output) 


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
    # DEBUG ----------------------------------------------------------------------------
    remix_units(
        path_to_nuts="build/data/administrative-borders-nuts.gpkg",
        path_to_gadm="build/data/administrative-borders-gadm.gpkg",
        path_to_output="build/data/national/units.geojson",
        resolution="national",
        all_countries=["Ireland", "United Kingdom"],
        layer_config={
            "Ireland": "nuts0",
            "United Kingdom": "nuts2",
            "Austria": "nuts0",
            "Belgium": "nuts0",
            "Bulgaria": "nuts0",
            "Croatia": "nuts0",
            "Cyprus": "nuts0",
            "Czech Republic": "nuts0",
            "Denmark": "nuts0",
            "Estonia": "nuts0",
            "Finland": "nuts0",
            "France": "nuts0",
            "Germany": "nuts0",
            "Greece": "nuts0",
            "Hungary": "nuts0",
            "Italy": "nuts0",
            "Latvia": "nuts0",
            "Lithuania": "nuts0",
            "Luxembourg": "nuts0",
            "Netherlands": "nuts0",
            "Poland": "nuts0",
            "Portugal": "nuts0",
            "Romania": "nuts0",
            "Slovakia": "nuts0",
            "Slovenia": "nuts0",
            "Spain": "nuts0",
            "Sweden": "nuts0",
            "Albania": "gadm0",
            "Bosnia and Herzegovina": "gadm0",
            "Macedonia, Republic of": "gadm0",
            "Montenegro": "gadm0",
            "Norway": "nuts0",
            "Serbia": "gadm0",
            "Switzerland": "nuts0",
        },
    )
    # ----------------------------------------------------------------------------------

    resolution = snakemake.wildcards[0]
    remix_units(
        path_to_nuts=snakemake.input.nuts,
        path_to_gadm=snakemake.input.gadm,
        path_to_output=snakemake.output[0],
        resolution=resolution,
        all_countries=snakemake.params.all_countries,
        layer_config=snakemake.params.layer_configs[resolution],
    )
