"""Take national potentials from JRC report and allocate to regions based on proxies."""
from enum import Enum
from pathlib import Path

import pandas as pd
import geopandas as gpd
import pycountry

PJ_TO_MWH = 1 / 3600 * 1e9
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


def biofuel_potential(paths_to_national_potentials, path_to_units, path_to_land_cover,
                      path_to_population, scenario, year, path_to_output):
    """Take national potentials from JRC report and allocate to regions based on proxies."""
    paths_to_national_potentials = [Path(path) for path in paths_to_national_potentials]
    national_potentials = [pd.read_csv(path, index_col=0, header=[0, 1])
                             .rename(index=eu_country_code_to_iso3)
                             .loc[:, (scenario, year)]
                             .rename(path.stem) * PJ_TO_MWH
                           for path in paths_to_national_potentials]
    units = gpd.read_file(path_to_units).set_index("id")
    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        national_potentials = [pd.Series(index=["EUR"], data=potential.sum(axis=0)).rename(potential.name)
                               for potential in national_potentials]
    total_potential = allocate_potentials(
        national_potentials=national_potentials,
        units=units,
        population=pd.read_csv(path_to_population, index_col=0)["population_sum"],
        land_cover=pd.read_csv(path_to_land_cover, index_col=0)
    )
    total_potential.to_csv(path_to_output, index=True, header=True)


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
    return sum(
        yessa(potential, units)
        for potential in national_potentials
    ).rename(NAME)


def yessa(national_potential, units):
    merged = units.merge(national_potential, right_index=True, left_on="country_code")
    return merged[PROXIES[national_potential.name]] * merged[national_potential.name]


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.

    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(eu_country_code) == 2, "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)
    if eu_country_code.lower() == "el":
        iso2 = "gr"
    elif eu_country_code.lower() == "uk":
        iso2 = "gb"
    elif eu_country_code.lower() == "bh": # this is just blantly wrong
        iso2 = "ba"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3


if __name__ == "__main__":
    biofuel_potential(
        paths_to_national_potentials=snakemake.input.national_potentials,
        path_to_units=snakemake.input.units,
        path_to_land_cover=snakemake.input.land_cover,
        path_to_population=snakemake.input.population,
        scenario=snakemake.params.scenario,
        year=snakemake.params.year,
        path_to_output=snakemake.output[0]
    )
