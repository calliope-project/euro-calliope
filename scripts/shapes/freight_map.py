import pandas as pd
import geopandas as gpd

from eurocalliopelib import utils

DRIVER = "GeoJSON"

FREIGHT_SECTORS = {
    "GT03": "FC_IND_MQ_E",
    "GT04": "FC_IND_FBT_E",
    "GT05": "FC_IND_TL_E",
    "GT06": "FC_IND_WP_E",
    "GT11": "FC_IND_MAC_E",
    "GT12": "FC_IND_TE_E",
    "GT13": "FC_IND_NSP_E"
}


def freight_map_creator(path_to_eu_freight, path_to_shapes, nuts_level, crs, industry_subsector, path_to_output):
    eu_freight = utils.read_eurostat_tsv(path_to_eu_freight).stack().to_xarray()

    eu_freight = (
        eu_freight
        .rename({"time": "year", "geo": "id"})
        .sortby("year")
        .squeeze()
        .drop("unit")
    )

    eu_freight = utils.rename_and_groupby(
        eu_freight, FREIGHT_SECTORS, dim_name="nst07", new_dim_name="cat_name"
    )

    eu_subsector_freight = eu_freight.sel(cat_name=industry_subsector)

    shapes = gpd.read_file(path_to_shapes).set_index("id").to_crs(crs)
    high_res_shapes = shapes[shapes.LEVL_CODE == nuts_level]
    eu_subsector_freight_reindexed = eu_subsector_freight.reindex({"id": high_res_shapes.index})

    # ASSUME: gross added value in a region without data in a given year is most
    # accurately represented by the nearest available year of data
    eu_subsector_freight_filled = eu_subsector_freight_reindexed.interpolate_na(
        "year", method="nearest", fill_value="extrapolate"
    )

    for year in eu_subsector_freight_filled.year.values:
        high_res_shapes[str(year)] = eu_subsector_freight_filled.sel(year=year).to_series()

    country_code_mapping = utils.convert_valid_countries(high_res_shapes.CNTR_CODE.unique())

    data_to_save = (
        high_res_shapes
        .rename(columns={"CNTR_CODE": "country_code"})
        .replace({"country_code": country_code_mapping})
        .loc[:, ["country_code", "geometry", *eu_subsector_freight_filled.year.astype(str).values]]
    )

    data_to_save.to_file(path_to_output, DRIVER)


if __name__ == "__main__":
    freight_map_creator(
        path_to_eu_freight=snakemake.input.eu_freight,
        path_to_shapes=snakemake.input.shapes,
        nuts_level=snakemake.params.nuts_level,
        crs=snakemake.params.crs,
        industry_subsector=snakemake.wildcards.industry_subsector,
        path_to_output=snakemake.output[0]
    )
