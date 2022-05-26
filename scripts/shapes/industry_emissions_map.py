import pandas as pd
import geopandas as gpd
from shapely import wkt

from eurocalliopelib import utils

DRIVER = "GeoJSON"

HOTMAPS_SECTORS = {
    "Iron and steel": "FC_IND_IS_E",
    "Cement": "FC_IND_NMM_E",
    "Paper and printing": "FC_IND_PPP_E",
    "Non-ferrous metals": "FC_IND_NFM_E",
    "Chemical industry": "FC_IND_CPC_E",
    "Glass": "FC_IND_NMM_E",
    "Other non-classified": "FC_IND_NSP_E",
    "Non-metallic mineral products": "FC_IND_NMM_E"
}


def geospatialise_industy_site_data(path_to_industry_sites, industry_subsector, crs, path_to_output):
    hotmaps_industry_sites = pd.read_csv(path_to_industry_sites, delimiter=";")
    hotmaps_industry_subsector_sites = (
        hotmaps_industry_sites
        .replace({"Subsector": HOTMAPS_SECTORS})
        .where(lambda x: x.Subsector == industry_subsector)
        .dropna(how="all")
    )
    spatial_points = hotmaps_industry_subsector_sites.geom.str.split(";", expand=True)[1].dropna().apply(wkt.loads)
    hotmap_industry_subsector_points = gpd.GeoDataFrame(
        hotmaps_industry_subsector_sites
        .drop("geom", axis=1)
        .loc[spatial_points.index],
        geometry=spatial_points,
        crs="epsg:4326"
    )
    hotmap_industry_subsector_points = hotmap_industry_subsector_points.to_crs(crs)

    hotmap_industry_subsector_points["Country"] = (
        hotmap_industry_subsector_points.Country.map(get_hotmaps_country_code)
    )
    # ASSUME: ETS data supersedes EPRTR data
    hotmap_industry_subsector_points["emissions"] = (
        hotmap_industry_subsector_points
        .Emissions_ETS_2014
        .fillna(hotmap_industry_subsector_points.Emissions_EPRTR_2014)
    )
    hotmap_industry_subsector_points = (
        hotmap_industry_subsector_points
        .where(hotmap_industry_subsector_points.emissions > 0)
        .dropna(subset=["emissions"])
    )
    hotmap_industry_subsector_points.to_file(path_to_output, DRIVER)


def get_hotmaps_country_code(country):
    if country == "Netherland":
        return utils.convert_country_code("netherlands")
    else:
        return utils.convert_country_code(country)


if __name__ == "__main__":
    geospatialise_industy_site_data(
        path_to_industry_sites=snakemake.input.industry_sites,
        crs=snakemake.params.crs,
        industry_subsector=snakemake.wildcards.industry_subsector,
        path_to_output=snakemake.output[0]
    )
