import pandas as pd
import pycountry


def scale_to_regional_resolution(df, region_country_mapping, populations):
    """
    Create regional electricity demand for controlled charging.
    ASSUME all road transport is subnationally distributed in proportion to population.
    """
    df_population_share = (
        populations.loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )

    regional_df = (
        pd.DataFrame(
            index=df.index,
            data={
                id: df[country_code]
                for id, country_code in region_country_mapping.items()
            },
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace(".", "-"))
    )
    pd.testing.assert_series_equal(regional_df.sum(axis=1), df.sum(axis=1))
    return regional_df


def scale_to_national_resolution(df):
    df.columns.name = None
    return df


def scale_to_continental_resolution(df):
    return df.sum(axis=1).to_frame("EUR")


def convert_annual_distance_to_electricity_demand(
    path_to_controlled_annual_demand: str,
    power_scaling_factor: float,
    first_year: int,
    final_year: int,
    conversion_factors: dict[str, float],
    country_codes: list[str],
):
    """
    Convert annual distance driven demand to electricity demand for
    controlled charging accounting for conversion factors.
    """
    df_energy_demand = (
        pd.read_csv(path_to_controlled_annual_demand, index_col=[1, 2])
        .xs(slice(first_year, final_year), level="year", drop_level=False)
        .assign(value=lambda x: x["value"] * x["vehicle_type"].map(conversion_factors))
        .groupby(["country_code", "year"])
        .sum()
        .loc[country_codes]
        .mul(power_scaling_factor)
        .squeeze()
        .unstack("country_code")
    )

    return -df_energy_demand


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution

    path_to_controlled_annual_demand = snakemake.input.annual_controlled_demand
    power_scaling_factor = snakemake.params.power_scaling_factor
    first_year = snakemake.params.first_year
    final_year = snakemake.params.final_year
    conversion_factors = snakemake.params.conversion_factors
    path_to_output = snakemake.output[0]
    country_codes = (
        [pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries],
    )
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )
    populations = pd.read_csv(snakemake.input.populations, index_col=0)

    df = convert_annual_distance_to_electricity_demand(
        path_to_controlled_annual_demand,
        power_scaling_factor,
        first_year,
        final_year,
        conversion_factors,
        country_codes,
    )

    if resolution == "continental":
        df = scale_to_continental_resolution(df)
    elif resolution == "national":
        df = scale_to_national_resolution(df)
    elif resolution in ["regional", "ehighways"]:
        df = scale_to_regional_resolution(
            df, region_country_mapping=region_country_mapping, populations=populations
        )
    else:
        raise ValueError("Input resolution is not recognised")
    df.T.to_csv(path_to_output, index_label=["id"])
