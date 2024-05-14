import pandas as pd
import pycountry


def scale_to_resolution_and_create_file(
    df, region_country_mapping, populations, resolution, output_path
):
    if resolution == "national":
        df = df
    elif resolution == "continental":
        df = df.sum(axis=1).to_frame("EUR")
    elif resolution in ["regional", "ehighways"]:
        df = scale_national_to_regional(df, region_country_mapping, populations)
    else:
        raise ValueError(f"Resolution {resolution} is not supported")
    df.tz_localize(None).rename_axis("utc-timestamp").to_csv(output_path)


def scale_national_to_regional(df, region_country_mapping, populations):
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


def get_national_ev_profiles(
    ev_profiles_path: str,
    dataset_name: str,
    demand_range: dict[str, int],
    first_year: int,
    final_year: int,
    country_neighbour_dict: dict[str, list[str]],
    country_codes: list[str],
):
    df_timeseries = (
        pd.read_csv(ev_profiles_path, index_col=[0, 1, 2], parse_dates=[0])
        .xs(slice(first_year, final_year), level="year")
        .unstack("country_code")
        .droplevel(level=0, axis="columns")
    )
    if "demand" in dataset_name:
        # Normalise demand and create min-max-equals timeseries
        df = (
            df_timeseries.groupby(by=lambda idx: idx.year)
            .transform(lambda x: x / x.sum())
            .mul(demand_range[dataset_name.split("-")[-1]])
        )
    elif "plugin" in dataset_name:
        # plugin-profiles are already normalised
        df = df_timeseries
    return df.pipe(fill_empty_country, country_neighbour_dict).loc[:, country_codes]


def fill_empty_country(df, country_neighbour_dict):
    for country, neighbours in country_neighbour_dict.items():
        assert country not in df.columns
        df[country] = df[neighbours].mean(axis=1)
    return df


if __name__ == "__main__":
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )
    populations = pd.read_csv(snakemake.input.populations, index_col=0)

    df = get_national_ev_profiles(
        ev_profiles_path=snakemake.input.ev_profiles,
        dataset_name=snakemake.wildcards.dataset_name,
        demand_range=snakemake.params.demand_range,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        country_neighbour_dict=snakemake.params.country_neighbour_dict,
        country_codes=([
            pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
        ]),
    )

    scale_to_resolution_and_create_file(
        df=df,
        region_country_mapping=region_country_mapping,
        populations=populations,
        resolution=snakemake.wildcards.resolution,
        output_path=snakemake.output[0],
    )
