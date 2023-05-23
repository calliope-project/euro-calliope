"""Take national potentials, allocate to regions based on feedstocks and proxies, and sum over feedstocks."""
from enum import Enum

import pandas as pd

PJ_TO_MWH = 1 / 3600 * 1e9
GJ_TO_MWH = 1 / 3600 * 1e3
NAME = "biofuel_potential_mwh_per_year"


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
                      path_to_population, scenario, potential_year, cost_year, proxies, paths_to_output):
    """Take national potentials from JRC report and allocate to regions based on proxies."""
    national_potentials = (
        pd
        .read_csv(path_to_national_potentials, index_col=["year", "scenario", "country_code", "feedstock"])["value"]
        .mul(PJ_TO_MWH)
        .xs((potential_year, scenario), level=("year", "scenario"))
    )
    national_costs = (
        pd
        .read_csv(path_to_national_costs, index_col=["year", "scenario", "country_code", "feedstock"])["value"]
        .div(GJ_TO_MWH)
        .xs((potential_year, scenario), level=("year", "scenario"))
    )
    units = pd.read_csv(path_to_units).set_index("id")
    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        cost = (
            national_costs
            .mul(national_potentials)
            .div(national_potentials.sum())
            .sum()
        )
        national_costs = pd.Series([cost], index=["EUR"]).rename_axis(index="country_code")
        national_potentials = (
            national_potentials
            .rename(lambda x: "EUR", level="country_code")
            .groupby(level=["country_code", "feedstock"])
            .sum()
        )

    total_potential = allocate_potentials(
        national_potentials=national_potentials,
        units=units,
        population=pd.read_csv(path_to_population, index_col=0).reindex(index=units.index)["population_sum"],
        land_cover=pd.read_csv(path_to_land_cover, index_col=0).reindex(index=units.index),
        proxies=proxies
    )
    total_potential.to_csv(paths_to_output.potentials, index=True, header=True)
    weighted_cost = (
        national_costs
        .mul(national_potentials)
        .div(national_potentials.sum())
        .sum()
    )
    with open(paths_to_output.costs, "w") as f_cost:
        f_cost.write(str(weighted_cost))


def allocate_potentials(national_potentials, units, population, land_cover, proxies):
    regional_potentials = (
        pd.merge(
            national_potentials.reset_index(),
            units["country_code"].reset_index(),
            on="country_code"
        )
        .pivot(index="id", columns="feedstock", values="value")
    )
    shares = (
        pd.concat(
            [population, land_cover[FOREST].sum(axis=1), land_cover[FARM].sum(axis=1)],
            axis=1,
            keys=['population', 'forest', "farmland"]
        )
        .groupby(units.country_code)
        .transform(lambda x: x / x.sum())
    )
    regional_potentials_sum = regional_potentials.groupby(proxies, axis=1).sum()
    regional_potentials_weighted = regional_potentials_sum.mul(shares).sum(axis=1)
    return regional_potentials_weighted.rename(NAME)


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
        proxies=snakemake.params.proxies,
        paths_to_output=snakemake.output
    )
