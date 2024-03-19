import numpy as np
import pandas as pd
import util


def national_to_continental_resolution(
        path_to_annual_demand: str,
        path_to_locations: str,
        path_to_populations: str
) -> pd.DataFrame:

    
    breakpoint()
    #    .groupby(by=['year', 'unit', 'end_use', 'cat_name'])
    #    .sum()
    #    .rename({'0': 'EUR'}, axis=1)
    # )

    return demand


def national_to_national_resolution(
        path_to_annual_demand: str,
        path_to_locations: str,
        path_to_populations: str
) -> pd.DataFrame:
    # is this really necessary?
    
    annual_demand = pd.read_csv(path_to_annual_demand, index_col=[0, 1, 2, 3])
    breakpoint()
    demand = annual_demand

    return demand


def national_to_regional_resolution(
        path_to_annual_demand: str,
        path_to_locations: str,
        path_to_populations: str
) -> pd.DataFrame:
    breakpoint()

    annual_demand = pd.read_csv(path_to_annual_demand, index_col=[0, 1, 2, 3])
    locations = pd.read_csv(path_to_locations, index_col=0)
    populations = pd.read_csv(path_to_populations, index_col=0)

    demand = 0
    return demand


def get_population_intensity(path_to_population, units):
    population_df = (
        pd.read_csv(path_to_population, index_col=0)
        .set_index(units.set_index(['id', 'country_code']).index)
    )
    return (
        population_df.div(population_df.sum(level='country_code')).population_sum
    )


def align_and_scale(orig_df, scaling_df, units):

    aligned_df = scaling_df.align(orig_df)
    scaled_df = aligned_df[0].mul(aligned_df[1]).dropna()

    # make sure numbers add up
    check_scaling_df = scaled_df.sum(level=orig_df.index.names)

    assert np.allclose(
        check_scaling_df.sum(level='country_code'),
        orig_df.reindex(check_scaling_df.index).sum(level='country_code'),
        equal_nan=True
    )

    return scaled_df


def subnational_pop_weighted_demand(
    units, heat_demand, heat_electricity_consumption, population
):
    concat_dfs = []
    population_intensity = get_population_intensity(population, units)

    concat_dfs.append(align_and_scale(
        heat_demand.xs('household', level='cat_name'), population_intensity, units
    ))

    concat_dfs.append(align_and_scale(
        heat_electricity_consumption.xs('household', level='cat_name'),
        population_intensity, units
    ))

    return concat_dfs


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution
    path_to_annual_demand = snakemake.input.annual_demand
    path_to_locations = snakemake.params.locations
    path_to_populations = snakemake.params.populations
    path_to_demand = snakemake.output.demand

    if resolution == "continental":
        demand = national_to_continental_resolution(
            path_to_annual_demand,
            path_to_locations,
            path_to_populations,
        )
    elif resolution == "national":
        demand = national_to_national_resolution(
            path_to_annual_demand,
            path_to_locations,
            path_to_populations,
        )
    elif resolution == "regional":
        demand = national_to_regional_resolution(
            path_to_annual_demand,
            path_to_locations,
            path_to_populations,
        )
    else:
        raise ValueError(f"Unknown resolution {resolution}")

    demand.to_csv(path_to_demand)
