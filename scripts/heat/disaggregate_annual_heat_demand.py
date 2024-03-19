import numpy as np
import pandas as pd
import util


def national_to_continental_resolution(
        path_to_annual_demand: str,
) -> pd.DataFrame:

    demand = (
        fill_annual_heat_demand_missing_values(path_to_annual_demand)
        .sum(axis=1)
        .to_frame('EUR')
    )

    return demand


def national_to_national_resolution(
        path_to_annual_demand: str,
) -> pd.DataFrame:

    demand = (
        fill_annual_heat_demand_missing_values(path_to_annual_demand)
    )

    return demand


def national_to_regional_resolution(
        path_to_annual_demand: str,
        path_to_locations: str,
        path_to_populations: str,
) -> pd.DataFrame:

    # TODO regional allocation is according to population, maybe we want to have it per GVA for commercial demand

    demand = (
        fill_annual_heat_demand_missing_values(path_to_annual_demand)
    )

    region_country_mapping = (
        pd.read_csv(path_to_locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )

    df_population_share = (
        pd.read_csv(path_to_populations, index_col=0)
        .loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )

    regional_df = (  # FIXME here we remove GEO,ISL,MDA,MLT,UKR from the demand data if we use the mapping from the default config
        pd
        .DataFrame(
            index=demand.index,
            data={
                id: demand[country_code]
                for id, country_code in region_country_mapping.items()
            },
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace(".", "-"))
    )

    return regional_df


def fill_annual_heat_demand_missing_values(  # FIXME we shouldn't fill NaNs with 0 because this means 0 demand, this is a temporary fix to avoid missing data issues
        path_to_annual_demand: str,
) -> pd.DataFrame:
    demand = (
        pd
        .read_csv(path_to_annual_demand, index_col=[0, 1, 2, 3, 4])
        .unstack("country_code")
        .pipe(lambda x: x.set_axis(x.columns.droplevel(level=0), axis=1))
        .fillna(0)  # this is where the fix needs to happen
    )
    return demand


def get_population_intensity(region_country_mapping, path_to_populations):
    region_country_mapping = (
        pd.read_csv(region_country_mapping, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )
    df_population_share = (
        pd.read_csv(path_to_populations, index_col=0)
        .loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )
    return df_population_share


def get_population_intensity_(path_to_population, units):
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

if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution
    path_to_annual_demand = snakemake.input.annual_demand
    path_to_locations = snakemake.params.locations
    path_to_populations = snakemake.params.populations
    path_to_demand = snakemake.output.demand

    if resolution == "continental":
        demand = national_to_continental_resolution(
            path_to_annual_demand,
        )
    elif resolution == "national":
        demand = national_to_national_resolution(
            path_to_annual_demand,
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
