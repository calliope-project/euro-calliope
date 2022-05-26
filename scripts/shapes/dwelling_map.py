
import geopandas as gpd

from eurocalliopelib import utils

DRIVER = "GeoJSON"


def dwelling_map_creator(path_to_dwellings, path_to_shapes, nuts_level, crs, path_to_output):
    dwellings = utils.read_eurostat_tsv(path_to_dwellings).stack().to_xarray()
    dwellings = dwellings.rename({"geo": "id"}).squeeze().drop(['unit', 'time'])

    shapes = gpd.read_file(path_to_shapes).set_index("id").to_crs(crs)
    high_res_shapes = shapes[shapes.LEVL_CODE == nuts_level]

    dwelling_density = dwellings.reindex({"id": high_res_shapes.index}) / high_res_shapes.area.to_xarray()

    # DW == all dwellings, "RES1" = single residence, "RES2" = 2 family residence, "RES_GE3" = residence with >= families
    # ASSUME: RES2 a multifamily home, not a single family home
    high_res_shapes["singlefamily_dwelling_density"] = dwelling_density.sel(housing="DW", building="RES1").to_series()
    high_res_shapes["multifamily_dwelling_density"] = (
        dwelling_density
        .sel(housing="DW", building=["RES2", "RES_GE3"])
        .sum("building", min_count=1)
        .to_series()
    )

    country_code_mapping = utils.convert_valid_countries(high_res_shapes.CNTR_CODE.unique())

    data_to_save = (
        high_res_shapes
        .rename(columns={"CNTR_CODE": "country_code"})
        .replace({"country_code": country_code_mapping})
        .loc[:, ["singlefamily_dwelling_density", "multifamily_dwelling_density", "country_code", "geometry"]]
    )
    data_to_save.to_file(path_to_output, DRIVER)


if __name__ == "__main__":
    dwelling_map_creator(
        path_to_dwellings=snakemake.input.dwellings,
        path_to_shapes=snakemake.input.shapes,
        nuts_level=snakemake.params.nuts_level,
        crs=snakemake.params.crs,
        path_to_output=snakemake.output[0]
    )
