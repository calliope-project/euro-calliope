import pandas as pd
import pycountry


def create_continental_timeseries(paths_to_input: list[str], country_codes: list[str]) -> pd.DataFrame:
    ts = create_national_timeseries(paths_to_input, country_codes)
    return ts.sum(axis=1).rename("EUR")


def create_national_timeseries(paths_to_input: list[str], country_codes: list[str]) -> pd.DataFrame:
    all_ts = [
        pd.read_csv(path, index_col='utc-timestamp', parse_dates=True)
        for path in paths_to_input
    ]
    return sum(all_ts).loc[:, country_codes]


def create_regional_timeseries(
    paths_to_input: list[str],
    country_codes: list[str],
    country_region_mapping: str,
    population: str,
) -> pd.DataFrame:
    """Create regional timeseries by 
        1. disaggregating the national timeseries into regional
        2. scaling each region according to the relative population in that region

    Output: a dataframe with columns over regions and rows over timestamps.
    """

    df_national = create_national_timeseries(paths_to_input, country_codes)

    country_region_mapping = (
        pd.read_csv(country_region_mapping, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )

    df_population_share = (
        pd.read_csv(population, index_col=0)
        .loc[:, "population_sum"]
        .reindex(country_region_mapping.keys())
        .groupby(by=country_region_mapping)
        .transform(lambda df: df / df.sum())
    )

    df_regional = (
        pd.DataFrame(
            index=df_national.index,
            data={
                id: df_national[country_code]
                for id, country_code in country_region_mapping.items()
            }
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace('.', '-'))    
    )

    pd.testing.assert_series_equal(df_regional.sum(axis=1), df_national.sum(axis=1))

    return df_regional


if __name__ == "__main__":
    resolution=snakemake.wildcards.resolution
    paths_to_input=snakemake.input.time_series
    country_codes=[pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries]
    path_to_output=snakemake.output[0]
    path_to_locations=snakemake.input.locations
    path_to_populations=snakemake.input.populations

    if resolution == "continental":
        ts = create_continental_timeseries(paths_to_input, country_codes)
    elif resolution == "national":
        ts = create_national_timeseries(paths_to_input, country_codes)
    elif resolution == "regional":
        ts = create_regional_timeseries(paths_to_input, country_codes, path_to_locations, path_to_populations)
    else:
        raise ValueError("Input resolution not recognised.")

    ts.to_csv(path_to_output)

