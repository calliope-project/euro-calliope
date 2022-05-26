import geopandas as gpd

from eurocalliopelib import utils

DRIVER = "GeoJSON"

EMPLOYEE_SECTORS = {  # Manual allocation to Eurostat subsectors, where some do not overlap well.
    "B": "FC_IND_MQ_E",  # Mining and quarrying
    "C10": "FC_IND_FBT_E",  # Manufacture of food products
    "C11": "FC_IND_FBT_E",  # Manufacture of beverages
    "C12": "FC_IND_FBT_E",  # Manufacture of tobacco products
    "C13": "FC_IND_TL_E",  # Manufacture of textiles
    "C14": "FC_IND_TL_E",  # Manufacture of wearing apparel
    "C15": "FC_IND_TL_E",  # Manufacture of leather and related products
    "C16": "FC_IND_WP_E",  # Manufacture of wood and of products of wood and cork, except furniture; manufacture of articles of straw and plaiting materials
    "C17": "FC_IND_PPP_E",  # Manufacture of paper and paper products
    "C18": "FC_IND_PPP_E",  # Printing and reproduction of recorded media
    "C20": "FC_IND_CPC_E",  # Manufacture of chemicals and chemical products
    "C21": "FC_IND_CPC_E",  # Manufacture of basic pharmaceutical products and pharmaceutical preparations
    "C22": "FC_IND_NSP_E",  # Manufacture of rubber and plastic products
    "C23": "FC_IND_NMM_E",  # Manufacture of other non-metallic mineral products
    "C24": "FC_IND_IS_E",  # Manufacture of basic metals  ASSUME: used for Iron and Steel
    "C25": "FC_IND_NFM_E",  # Manufacture of fabricated metal products, except machinery and equipment  ASSUME: used for non-iron and steel
    "C26": "FC_IND_MAC_E",  # Manufacture of computer, electronic and optical products
    "C27": "FC_IND_MAC_E",  # Manufacture of electrical equipment
    "C28": "FC_IND_MAC_E",  # Manufacture of machinery and equipment n.e.c.
    "C29": "FC_IND_TE_E",  # Manufacture of motor vehicles, trailers and semi-trailers
    "C30": "FC_IND_TE_E",  # Manufacture of other FC_IND_TE_E
    "C31": "FC_IND_NSP_E",  # Manufacture of furniture
    "C32": "FC_IND_NSP_E",  # Other manufacturing
    "F": "FC_IND_CON_E",  # Construction
}


def employees_map_creator(path_to_eu_business_statistics, path_to_shapes, nuts_level, crs, industry_subsector, path_to_output):
    eu_business_statistics = utils.read_eurostat_tsv(path_to_eu_business_statistics).stack().to_xarray()

    eu_business_statistics = eu_business_statistics.rename({"time": "year", "geo": "id"}).sortby("year")

    # Options to select on: V11210 == local units,  V16110 == persons employed, V13320 == wages and salaries
    # ASSUME: persons employed is the best indicator of industry activity

    eu_employees = eu_business_statistics.sel(indic_sb="V16110")
    eu_employees = utils.rename_and_groupby(
        eu_employees, EMPLOYEE_SECTORS, dim_name="nace_r2", new_dim_name="cat_name")

    eu_subsector_employees = eu_employees.sel(cat_name=industry_subsector)

    shapes = gpd.read_file(path_to_shapes).set_index("id").to_crs(crs)
    high_res_shapes = shapes[shapes.LEVL_CODE == nuts_level]
    eu_subsector_employees_reindexed = eu_subsector_employees.reindex({"id": high_res_shapes.index})

    # ASSUME: gross added value in a region without data in a given year is most
    # accurately represented by the nearest available year of data
    eu_subsector_employees_filled = eu_subsector_employees_reindexed.interpolate_na(
        "year", method="nearest", fill_value="extrapolate"
    )

    for year in eu_subsector_employees_filled.year.values:
        high_res_shapes[str(year)] = eu_subsector_employees_filled.sel(year=year).to_series()

    country_code_mapping = utils.convert_valid_countries(high_res_shapes.CNTR_CODE.unique())

    data_to_save = (
        high_res_shapes
        .rename(columns={"CNTR_CODE": "country_code"})
        .replace({"country_code": country_code_mapping})
        .loc[:, ["country_code", "geometry", *eu_subsector_employees_filled.year.astype(str).values]]
    )

    data_to_save.to_file(path_to_output, DRIVER)


if __name__ == "__main__":
    employees_map_creator(
        path_to_eu_business_statistics=snakemake.input.eu_business_statistics,
        path_to_shapes=snakemake.input.shapes,
        nuts_level=snakemake.params.nuts_level,
        crs=snakemake.params.crs,
        industry_subsector=snakemake.wildcards.industry_subsector,
        path_to_output=snakemake.output[0]
    )
