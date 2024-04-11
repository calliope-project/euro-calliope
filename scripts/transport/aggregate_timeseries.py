import pandas as pd


def create_continental_timeseries(paths_to_input: list[str]) -> pd.DataFrame:
    ts = create_national_timeseries(paths_to_input)
    return ts.sum(axis=1).rename("EUR")


def create_national_timeseries(paths_to_input: list[str]) -> pd.DataFrame:
    all_ts = [
        pd.read_csv(path, index_col="utc-timestamp", parse_dates=True)
        for path in paths_to_input
    ]
    return sum(all_ts)


def create_regional_timeseries(
    paths_to_input: list[str],
    region_country_mapping: str,
    population: str,
) -> pd.DataFrame:
    """Create regional timeseries by
        1. disaggregating the national timeseries into regional
        2. scaling each region according to the relative population in that region

    Output: a dataframe with columns over regions and rows over timestamps.

    ASSUME all road transport is subnationally distributed in proportion to population.
    """

    df_national = create_national_timeseries(paths_to_input)

    region_country_mapping = (
        pd.read_csv(region_country_mapping, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )

    df_population_share = (
        pd.read_csv(population, index_col=0)
        .loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )

    df_regional = (
        pd.DataFrame(
            index=df_national.index,
            data={
                id: df_national[country_code]
                for id, country_code in region_country_mapping.items()
            },
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace(".", "-"))
    )

    pd.testing.assert_series_equal(df_regional.sum(axis=1), df_national.sum(axis=1))

    return df_regional


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution
    paths_to_input = snakemake.input.time_series
    path_to_output = snakemake.output[0]
    path_to_locations = snakemake.input.locations
    path_to_populations = snakemake.input.populations

    if resolution == "continental":
        ts = create_continental_timeseries(paths_to_input)
    elif resolution == "national":
        ts = create_national_timeseries(paths_to_input)
    elif resolution == "regional":
        ts = create_regional_timeseries(
            paths_to_input, path_to_locations, path_to_populations
        )
    elif resolution == "ehighways":
        ts = create_regional_timeseries(
            paths_to_input, path_to_locations, path_to_populations
        )
    else:
        raise ValueError("Input resolution not recognised.")

    ts.to_csv(path_to_output)
