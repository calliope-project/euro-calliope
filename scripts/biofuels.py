"""Take national potentials from JRC report and allocate to regions based on proxies."""
from enum import Enum
from pathlib import Path

import pandas as pd
import geopandas as gpd

from eurocalliopelib import utils

PJ_TO_MWH = 1 / 3600 * 1e9
GJ_TO_MWH = 1 / 3600 * 1e3
NAME = "biofuel_potential_mwh_per_year"

SCENARIOS = {
    "low": "Low availability scenario",
    "medium": "Medium availability scenario",
    "high": "High availability scenario"
}

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


def biofuel_potential(paths_to_national_potentials, paths_to_costs, path_to_units, path_to_land_cover,
                      path_to_population, scenario, potential_year, cost_year, path_to_potentials, path_to_costs):
    """Take national potentials from JRC report and allocate to regions based on proxies."""
    scenario = SCENARIOS[scenario]
    paths_to_national_potentials = [Path(path) for path in paths_to_national_potentials]
    paths_to_costs = [Path(path) for path in paths_to_costs]
    national_potentials = pd.concat(
        [pd.read_csv(path, index_col=0, header=[0, 1])
           .rename(index=utils.eu_country_code_to_iso3)
           .loc[:, (scenario, potential_year)]
           .rename(path.stem) * PJ_TO_MWH
         for path in paths_to_national_potentials],
        axis=1
    )
    costs = pd.concat(
        [pd.read_csv(path, index_col=0)
           .rename(index=utils.eu_country_code_to_iso3)
           .loc[:, cost_year]
           .rename(path.stem) / GJ_TO_MWH
         for path in paths_to_costs],
        axis=1
    )
    if scenario == "Low availability scenario":
        costs = costs * 1.1
    elif scenario == "High availability scenario":
        costs = costs * 0.9
    units = gpd.read_file(path_to_units).set_index("id")
    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        national_potentials = pd.DataFrame(index=["EUR"], data=national_potentials.sum(axis=0).to_dict())
        costs = pd.DataFrame(index=["EUR"], data=costs.mean(axis=0).to_dict())

    total_potential = allocate_potentials(
        national_potentials=national_potentials,
        units=units,
        population=pd.read_csv(path_to_population, index_col=0).reindex(index=units.index)["population_sum"],
        land_cover=pd.read_csv(path_to_land_cover, index_col=0).reindex(index=units.index)
    )
    total_potential.to_csv(path_to_potentials, index=True, header=True)
    weighted_cost = (costs * national_potentials / national_potentials.sum().sum()).sum().sum()
    with open(path_to_costs, "w") as f_cost:
        f_cost.write(str(weighted_cost))


def allocate_potentials(national_potentials, units, population, land_cover):
    units = pd.concat([
        units,
        population.groupby(units.country_code)
                  .transform(lambda x: x / x.sum())
                  .rename("population_share"),
        land_cover[FOREST].sum(axis=1)
                          .groupby(units.country_code)
                          .transform(lambda x: x / x.sum())
                          .rename("forest_share"),
        land_cover[FARM].sum(axis=1)
                        .groupby(units.country_code)
                        .transform(lambda x: x / x.sum())
                        .rename("farm_share"),
    ], axis=1)
    units = units.merge(national_potentials, right_index=True, left_on="country_code")
    return sum(
        units[PROXIES[potential]] * units[potential]
        for potential in national_potentials
    ).rename(NAME)


if __name__ == "__main__":
    biofuel_potential(
        paths_to_national_potentials=snakemake.input.national_potentials,
        paths_to_costs=snakemake.input.costs,
        path_to_units=snakemake.input.units,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_population=snakemake.input.population,
        scenario=snakemake.wildcards.scenario,
        potential_year=snakemake.params.potential_year,
        cost_year=snakemake.params.cost_year,
        path_to_potentials=snakemake.output.potentials,
        path_to_costs=snakemake.output.costs
    )
