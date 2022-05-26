import pandas as pd
import geopandas as gpd

from eurocalliopelib import utils

DRIVER = "GeoJSON"


def gva_map_creator(path_to_eu_gva, path_to_ch_gva, path_to_shapes, nuts_level, crs, path_to_output):
    eu_gva = utils.read_eurostat_tsv(path_to_eu_gva).stack().to_xarray()

    # ASSUME: economic activity is a better indicator than
    eu_gva = eu_gva.sel(currency="MIO_EUR")

    shapes = gpd.read_file(path_to_shapes).set_index("id").to_crs(crs)
    high_res_shapes = shapes[shapes.LEVL_CODE == nuts_level]
    eu_gva = eu_gva.rename({"time": "year", "geo": "id"}).sortby("year")

    eu_gva_reindexed = eu_gva.reindex({"id": high_res_shapes.index})

    # KNOWN-ISSUE: We use NUTS2016 shapes, but Norway's data is based on older regions,
    # which means that NO060 (a reasonably large area) is lost without this manual gap filling
    eu_gva_reindexed.loc[{"id": "NO060"}] = eu_gva.sel(id=["NO061", "NO062"]).sum("id", min_count=1)

    # ASSUME: NACE Revision 2 commercial subsectors are ["G-J", "K-N", "O-U"]
    eu_gva_commercial = eu_gva_reindexed.sel(nace_r2=["G-J", "K-N", "O-U"]).sum("nace_r2", min_count=1)

    eu_gva_total = eu_gva_reindexed.sel(nace_r2="TOTAL")

    # ASSUME: total GVA can be used as a placeholder for commercial sector GVA when
    # commercial sector GVA data is not available
    eu_gva_filled = eu_gva_commercial.fillna(eu_gva_total)

    ch_gva = pd.read_csv(path_to_ch_gva, index_col=[0, 1], squeeze=True).to_xarray()

    eu_gva_filled = eu_gva_filled.fillna(ch_gva)
    assert eu_gva_filled.sel(id=ch_gva.id).notnull().any()

    # ASSUME: gross added value in a region without data in a given year is most
    # accurately represented by the nearest available year of data
    eu_gva_filled = eu_gva_filled.interpolate_na("year", method="nearest", fill_value="extrapolate")

    for year in eu_gva_filled.year.values:
        high_res_shapes[str(year)] = eu_gva_filled.sel(year=year).to_series()

    country_code_mapping = utils.convert_valid_countries(high_res_shapes.CNTR_CODE.unique())

    data_to_save = (
        high_res_shapes
        .rename(columns={"CNTR_CODE": "country_code"})
        .replace({"country_code": country_code_mapping})
        .loc[:, ["country_code", "geometry", *eu_gva_filled.year.astype(str).values]]
    )

    data_to_save.to_file(path_to_output, DRIVER)


if __name__ == "__main__":
    gva_map_creator(
        path_to_eu_gva=snakemake.input.eu_gva,
        path_to_ch_gva=snakemake.input.ch_gva,
        path_to_shapes=snakemake.input.shapes,
        nuts_level=snakemake.params.nuts_level,
        crs=snakemake.params.crs,
        path_to_output=snakemake.output[0]
    )
