"""Remixes NUTS, LAU, and GADM data to form the units of the analysis."""

import fiona
import geopandas as gpd
import pandas as pd
from eurocalliopelib import utils

DRIVER = "GeoJSON"


def remix_units(
    path_to_nuts: str,
    path_to_gadm: str,
    path_to_ehighways: str,
    path_to_output: str,
    resolution: str,
    resolution_config: dict,
    all_countries: list[str],
):
    """Remixes NUTS, LAU, and GADM data to form the units of the analysis."""
    source_layers = _read_source_layers(path_to_nuts, path_to_gadm, path_to_ehighways)
    _validate_source_layer_crs(source_layers)
    _validate_resolution_config(all_countries, resolution_config, resolution)
    units = _build_layer(resolution_config, source_layers)
    _verify_all_countries_captured(all_countries, units, resolution)
    if resolution == "continental":
        # sum all of the units together into one block
        units = _merge_national_shapes_to_continental_layer(units)
    elif resolution == "ehighways":
        units = _rename_ehighways_countries(units)
    units.to_file(path_to_output, driver=DRIVER)


def _build_layer(
    resolution_config: dict[str, str], source_layers: dict[str, gpd.GeoDataFrame]
):
    """
    Inputs:
    - resolution_config: a mapping of which statistical unit resolution each country should use as its base
    - source_layers: a dict with keys as each statistical unit type (nuts[0-3] or gadmn[0-3])
        and values as a GeoDataFrame representing the collection of those units
    Returns:
    - layer: the GeoDataFrame covering the whole scope in source_layers, where each
        country is formed from units corresponding to the resolution specified in
        resolution_config
    """
    crs = [layer.crs for layer in source_layers.values()][0]
    layer = pd.concat([
        source_layers[source_layer][
            source_layers[source_layer].country_code == _iso3(country)
        ]
        for country, source_layer in resolution_config.items()
    ])
    assert isinstance(layer, pd.DataFrame)
    return gpd.GeoDataFrame(layer, crs=crs).reset_index(drop=True)


def _read_source_layers(
    path_to_nuts: str, path_to_gadm: str, path_to_ehighways: str
) -> dict[str, gpd.GeoDataFrame]:
    """Returns:
    - source_layers: a dict with keys as each statistical unit type (nuts[0-3], gadmn[0-3], or ehighways)
        and values as a geodataframe representing the collection of those units.
    """

    source_layers = {
        layer_name: gpd.read_file(src, layer=layer_name)
        for src in [path_to_nuts, path_to_gadm, path_to_ehighways]
        for layer_name in fiona.listlayers(src)
    }
    return source_layers


def _validate_source_layer_crs(source_layers: dict[str, gpd.GeoDataFrame]):
    crs = [layer.crs for layer in source_layers.values()]
    assert not crs or crs.count(crs[0]) == len(
        crs
    ), "Source layers have different crs. They must match."


def _validate_resolution_config(
    all_countries: list[str], resolution_config: dict[str, str], resolution: str
):
    missing_countries = set(all_countries).difference(resolution_config)

    assert not missing_countries, f"Missing countries in {resolution} resolution configuration: {missing_countries}"


def _verify_all_countries_captured(
    countries: list[str], units: gpd.GeoDataFrame, resolution: str
):
    missing_countries = set(units.country_code).difference([
        _iso3(country) for country in countries
    ])
    assert (
        not missing_countries
    ), f"Countries are missing in {resolution} resolution shapes: {missing_countries}"


def _iso3(country_name: str) -> str:
    return utils.convert_country_code(country_name, output="alpha3")


def _merge_national_shapes_to_continental_layer(
    units: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    "For the continental resolution, dissolve all"
    units = units.assign(all=1).dissolve("all")
    index = units.index[0]
    units.loc[index, "id"] = "EUR"
    units.loc[index, "country_code"] = "EUR"
    units.loc[index, "name"] = "Europe"
    units.loc[index, "type"] = "continent"
    units.loc[index, "proper"] = 1
    return units


def _rename_ehighways_countries(units: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """Align non-ehighways IDs to match ehighways ones by appending "_1" to country codes, e.g. "BEL" -> "BEL_1".

    Method is future-proofed to handle higher resolution than country code (e.g. using GADM3 units in one country).
    It does this by ordering all non-ehighways IDs alphabetically and then turning them into "XXX_#" style codes.
    """
    order = (
        units.id.str.extractall(r"(\d+)")
        .groupby(level=0)
        .sum()
        .reindex(units.index)
        .fillna(1)
        .astype(int)
    )
    sorted_units = units.assign(order=order).sort_values(by=["country_code", "order"])

    new_names = (
        sorted_units.country_code
        + "_"
        + sorted_units.groupby("country_code").cumcount().add(1).astype(str)
    )

    sorted_units["id"] = sorted_units.id.where(
        sorted_units.type == "ehighways", new_names
    )

    return sorted_units.drop("order", axis=1)


if __name__ == "__main__":
    remix_units(
        path_to_nuts=snakemake.input.nuts,
        path_to_gadm=snakemake.input.gadm,
        path_to_ehighways=snakemake.input.ehighways,
        path_to_output=snakemake.output[0],
        all_countries=snakemake.params.all_countries,
        resolution_config=snakemake.params.resolution_config,
        resolution=snakemake.wildcards[0],
    )
