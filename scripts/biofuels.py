"""Take national potentials, allocate to regions based on feedstocks and proxies, and sum over feedstocks."""
from enum import Enum

import pandas as pd
import xarray as xr

PJ_TO_MWH = 1 / 3600 * 1e9
GJ_TO_MWH = 1 / 3600 * 1e3
NAME = "biofuel_potential_mwh_per_year"

PROXIES = {
    "forestry-energy-residues": "forest_share",
    "forestry-care-residues": "forest_share",
    "roundwood-chips": "forest_share",
    "roundwood-fuelwood": "forest_share",
    "secondary-forestry-residues-sawdust": "forest_share",
    "secondary-forestry-residues-woodchips": "forest_share",
    "manure": "farm_share",
    "primary-agricultural-residues": "farm_share",
    "municipal-waste": "population_share",
    "landscape-care-residues": "population_share",
    "sludge": "population_share"
}


class GlobCover(Enum):
    """Original categories taken from GlobCover 2009 land cover."""
    POST_FLOODING = "lc_11"
    RAINFED_CROPLANDS = "lc_14"
    MOSAIC_CROPLAND = "lc_20"
    MOSAIC_VEGETATION = "lc_30"
    CLOSED_TO_OPEN_BROADLEAVED_FOREST = "lc_40"
    CLOSED_BROADLEAVED_FOREST = "lc_50"
    OPEN_BROADLEAVED_FOREST = "lc_60"
    CLOSED_NEEDLELEAVED_FOREST = "lc_70"
    OPEN_NEEDLELEAVED_FOREST = "lc_90"
    CLOSED_TO_OPEN_MIXED_FOREST = "lc_100"
    MOSAIC_FOREST = "lc_110"
    MOSAIC_GRASSLAND = "lc_120"
    CLOSED_TO_OPEN_SHRUBLAND = "lc_130"
    CLOSED_TO_OPEN_HERBS = "lc_140"
    SPARSE_VEGETATION = "lc_150"
    CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST = "lc_160" # doesn't exist in Europe
    CLOSED_REGULARLY_FLOODED_FOREST = "lc_170" # doesn't exist in Europe
    CLOSED_TO_OPEN_REGULARLY_FLOODED_GRASSLAND = "lc_180" # roughly 2.3% of land in Europe
    ARTIFICAL_SURFACES_AND_URBAN_AREAS = "lc_190"
    BARE_AREAS = "lc_200"
    WATER_BODIES = "lc_210"
    PERMANENT_SNOW = "lc_220"
    NO_DATA = "lc_230"


FARM = [GlobCover.POST_FLOODING.value, GlobCover.RAINFED_CROPLANDS.value,
        GlobCover.MOSAIC_CROPLAND.value, GlobCover.MOSAIC_VEGETATION.value]
FOREST = [GlobCover.CLOSED_TO_OPEN_BROADLEAVED_FOREST.value, GlobCover.CLOSED_BROADLEAVED_FOREST.value,
          GlobCover.OPEN_BROADLEAVED_FOREST.value, GlobCover.CLOSED_NEEDLELEAVED_FOREST.value,
          GlobCover.OPEN_NEEDLELEAVED_FOREST.value, GlobCover.CLOSED_TO_OPEN_MIXED_FOREST.value,
          GlobCover.MOSAIC_FOREST.value, GlobCover.CLOSED_TO_OPEN_REGULARLY_FLOODED_FOREST.value,
          GlobCover.CLOSED_REGULARLY_FLOODED_FOREST.value]


def biofuel_potential(path_to_national_potentials, path_to_national_costs, path_to_units, path_to_land_cover,
                      path_to_population, scenario, potential_year, cost_year, paths_to_output):
    """Take national potentials from JRC report and allocate to regions based on proxies."""
    national_potentials = (
        pd
        .read_csv(path_to_national_potentials, index_col=["year", "scenario", "country_code", "feedstock"])["value"]
        .mul(PJ_TO_MWH)
        .to_xarray()
        .sel(year=potential_year, scenario=scenario)
    )
    national_costs = (
        pd
        .read_csv(path_to_national_costs, index_col=["year", "scenario", "country_code", "feedstock"])["value"]
        .div(GJ_TO_MWH)
        .to_xarray()
        .sel(year=cost_year, scenario=scenario)
    )
    units = pd.read_csv(path_to_units).set_index("id")
    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        national_costs = (
            (national_costs * national_potentials / national_potentials.sum(["country_code", "feedstock"]))
            .sum(["country_code", "feedstock"])
            .expand_dims({"country_code": ["EUR"]})
        )
        national_potentials = national_potentials.sum("country_code").expand_dims({"country_code": ["EUR"]})

    total_potential = allocate_potentials(
        national_potentials=national_potentials,
        units=units,
        population=pd.read_csv(path_to_population, index_col=0).reindex(index=units.index)["population_sum"],
        land_cover=pd.read_csv(path_to_land_cover, index_col=0).reindex(index=units.index)
    )
    total_potential.to_csv(paths_to_output.potentials, index=True, header=True)
    weighted_cost = (
        (national_costs * national_potentials / national_potentials.sum(["country_code", "feedstock"]))
        .sum(["country_code", "feedstock"])
        .item()
    )
    with open(paths_to_output.costs, "w") as f_cost:
        f_cost.write(str(weighted_cost))


def allocate_potentials(national_potentials, units, population, land_cover):
    ds = national_potentials.copy()
    ds = (
        pd
        .merge(units["country_code"].reset_index(), ds.to_series().reset_index(), on="country_code")
        .set_index(["id", "feedstock"])
        .to_xarray()
    )
    ds["population_share"] = (
        population
        .groupby(units.country_code)
        .transform(lambda x: x / x.sum())
        .rename("population_share")
    )
    ds["forest_share"] = (
        land_cover[FOREST]
        .sum(axis=1)
        .groupby(units.country_code)
        .transform(lambda x: x / x.sum())
        .rename("forest_share")
    )
    ds["farm_share"] = (
        land_cover[FARM]
        .sum(axis=1)
        .groupby(units.country_code)
        .transform(lambda x: x / x.sum())
        .rename("farm_share")
    )
    ds["proxy"] = pd.Series(PROXIES).rename("proxy").rename_axis(index="feedstock").to_xarray()
    proxy_value = xr.ones_like(ds.value)
    for feedstock in ds.feedstock:
        for id in ds.id:
            proxy = ds.sel(id=id, feedstock=feedstock).proxy.item()
            proxy_value.loc[dict(id=id, feedstock=feedstock)] = ds.sel(id=id)[proxy].item()
    ds["value"] = ds.value * proxy_value
    return ds.value.sum("feedstock").rename(NAME).to_series()


if __name__ == "__main__":
    biofuel_potential(
        path_to_national_potentials=snakemake.input.national_potentials,
        path_to_national_costs=snakemake.input.costs,
        path_to_units=snakemake.input.units,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_population=snakemake.input.population,
        scenario=snakemake.wildcards.scenario,
        potential_year=snakemake.params.potential_year,
        cost_year=snakemake.params.cost_year,
        paths_to_output=snakemake.output
    )
