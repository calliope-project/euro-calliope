import geopandas as gpd
import pandas as pd
import shapely.geometry
from eurocalliopelib import utils

DRIVER = "GeoJSON"


def remix_units(
    path_to_nuts, path_to_gadm, path_to_output, layer_name, layer_config, all_countries
):
    """Remixes NUTS, LAU, and GADM data to form the units of the analysis.
    source_layers: a dict with keys as each geographical layer type (nutsX or gadmnX)
    and values as a geodataframe of the POLYGONS
    """
    source_layers = _read_source_layers(path_to_nuts, path_to_gadm)
    _validate_source_layers(source_layers)
    _validate_layer_config(all_countries, layer_config, layer_name)
    layer = _build_layer(layer_config, source_layers)
    _validate_layer(layer, layer_name, all_countries)
    if layer_name == "continental":  # treat special case
        layer = _continental_layer(layer)
    _write_layer(layer, path_to_output)


def _read_source_layers(path_to_nuts, path_to_gadm):
    source_layers = {
        layer_name: gpd.read_file(path_to_nuts, layer=layer_name)
        for layer_name in fiona.listlayers(path_to_nuts)
    }
    source_layers.update(
        {
            layer_name: gpd.read_file(path_to_gadm, layer=layer_name)
            for layer_name in fiona.listlayers(path_to_gadm)
        }
    )
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


def _build_layer(country_to_source_map, source_layers):
    crs = [layer.crs for layer in source_layers.values()][0]
    layer = pd.concat(
        [
            source_layers[source_layer][
                source_layers[source_layer].country_code == _iso3(country)
            ]
            for country, source_layer in country_to_source_map.items()
        ]
    )
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


def dissolve_nuts3(
    # path_to_disaggregated_units,
    path_to_output,
    nuts_year,
    path_to_statistical_to_custom_units,
):
    """Dissolve NUTS3 data to custom regions."""

    stat_to_custom_df = pd.read_csv(path_to_statistical_to_custom_units, header=0)
    # Any clusters which are whole countries will be whole countries at the subregional scale too
    nuts_year_string = f"NUTS3_{nuts_year}"
    locations = (
        stat_to_custom_df.dropna(subset=[nuts_year_string])
        .set_index(
            stat_to_custom_df.dropna(subset=[nuts_year_string])[nuts_year_string]
        )
        .id.append(
            stat_to_custom_df[stat_to_custom_df.Source != "NUTS3"]
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
    # DEBUG CONFIGURATION --------------------------------------------------------------
    dissolve_nuts3(
        # path_to_disaggregated_units=snakemake.input.disaggregated_units, # no-longer needed
        path_to_statistical_to_custom_units="config/shapes/statistical-to-ehighways-units.csv",
        nuts_year=2006,
        path_to_output="build/data/ehighways/units.geojson",
    )
    # ----------------------------------------------------------------------------------

    dissolve_nuts3(
        # path_to_disaggregated_units=snakemake.input.disaggregated_units, # no-longer needed
        path_to_statistical_to_custom_units=snakemake.input.statistical_to_custom_units,
        nuts_year=int(snakemake.params.nuts_year),
        path_to_output=snakemake.output[0],
    )

    remix_units(
        path_to_nuts=snakemake.input.nuts,
        path_to_gadm=snakemake.input.gadm,
        path_to_output=snakemake.output[0],
        layer_name=layer_name,
        all_countries=snakemake.params.all_countries,
        layer_config=snakemake.params.layer_configs[layer_name],
    )
