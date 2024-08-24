import pandas as pd


def national_to_continental_resolution(
    annual_demand: pd.DataFrame, **kwargs
) -> pd.DataFrame:
    demand = annual_demand.sum(axis=1).to_frame("EUR")

    return demand


def national_to_national_resolution(
    annual_demand: pd.DataFrame, **kwargs
) -> pd.DataFrame:
    return annual_demand


def national_to_regional_resolution(
    annual_demand: pd.DataFrame,
    region_country_mapping: dict,
    populations: pd.DataFrame,
) -> pd.DataFrame:
    # ASSUME national heat demand is spatially distributed by population.
    # TODO maybe we want to have it per GVA for commercial demand
    df_population_share = (
        populations.loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )
    regional_df = (
        pd.DataFrame(
            index=annual_demand.index,
            data={
                id: annual_demand[country_code]
                for id, country_code in region_country_mapping.items()
            },
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace(".", "_"))
    )

    pd.testing.assert_series_equal(regional_df.sum(axis=1), annual_demand.sum(axis=1))

    return regional_df


def read_data(path_to_file: str):
    return (
        pd.read_csv(path_to_file, index_col=[0, 1, 2, 3])
        .squeeze()
        .unstack("country_code")
    )


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution
    annual_demand = read_data(snakemake.input.annual_demand)
    electrified = read_data(snakemake.input.electricity)
    populations = pd.read_csv(snakemake.input.populations, index_col=0)
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )

    if resolution == "continental":
        rescale_function = national_to_continental_resolution
    elif resolution == "national":
        rescale_function = national_to_national_resolution
    elif resolution in ["regional", "ehighways"]:
        rescale_function = national_to_regional_resolution
    else:
        raise ValueError(f"Unknown resolution {resolution}")

    (
        rescale_function(
            annual_demand,
            populations=populations,
            region_country_mapping=region_country_mapping,
        ).to_csv(snakemake.output.total_demand)
    )
    (
        rescale_function(
            electrified,
            populations=populations,
            region_country_mapping=region_country_mapping,
        ).to_csv(snakemake.output.electricity)
    )
