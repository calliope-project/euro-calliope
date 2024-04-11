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


def extract_national_ev_numbers(
    path_to_ev_numbers: str,
    transport_scaling_factor: float,
    first_year: int,
    final_year: int,
    conversion_factors: dict[str, float],
    battery_sizes: dict[str, float],
    country_codes: list[str],
):
    # Extract number of EVs per vehicle type
    df_ev_numbers = (
        pd.read_csv(path_to_ev_numbers, index_col=[0, 1, 2, 3, 4])
        .squeeze()
        .droplevel(["vehicle_subtype", "section"])
    )
    assert df_ev_numbers.isnull().values.any() == False
    # Compute max. distance travelled per full battery for one EV [in Mio km / vehicle]
    battery_size = pd.DataFrame.from_dict(
        {
            vehicle: battery_sizes[vehicle] / conversion_factors[vehicle]
            for vehicle in battery_sizes
        },
        orient="index",
        columns=["value"],
    )

    # Compute available chargeable distance per vehicle type [in transport scaling unit km]
    df_ev_chargeable_distance = (
        df_ev_numbers.align(battery_size, level="vehicle_type")[1]
        .squeeze()
        .mul(df_ev_numbers)
        .groupby(level=["country_code", "year"])
        .sum()
        .loc[country_codes]
        .unstack("year")
        .mul(transport_scaling_factor)
    )
    df_ev_chargeable_distance = df_ev_chargeable_distance.assign(**{
        str(year): df_ev_chargeable_distance[2015] for year in range(2016, 2019)
    })

    df_ev_chargeable_distance.columns = df_ev_chargeable_distance.columns.astype(int)
    df_ev_chargeable_distance.index.name = None

    return df_ev_chargeable_distance[range(first_year, final_year + 1)].T


def compute_weighted_share_per_vehicle_type(path_to_vehicle_type_distance):
    return (
        pd.read_csv(path_to_vehicle_type_distance, index_col=[0, 1, 2, 3, 4])
        .droplevel(["vehicle_subtype", "section"])
        .groupby(level=["country_code", "year"])
        .transform(lambda x: (x / x.sum()))
    )


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution

    path_to_ev_numbers = snakemake.input.ev_vehicle_number
    transport_scaling_factor = snakemake.params.transport_scaling_factor
    first_year = snakemake.params.first_year
    final_year = snakemake.params.final_year
    conversion_factors = snakemake.params.conversion_factors
    battery_sizes = snakemake.params.battery_sizes
    path_to_output = snakemake.output[0]
    country_codes = [
        pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
    ]
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )
    populations = pd.read_csv(snakemake.input.populations, index_col=0)

    df = extract_national_ev_numbers(
        path_to_ev_numbers,
        transport_scaling_factor,
        first_year,
        final_year,
        conversion_factors,
        battery_sizes,
        country_codes,
    )
    breakpoint()
    if resolution == "continental":
        df = scale_to_continental_resolution(df)
    elif resolution == "national":
        df = scale_to_national_resolution(df)
    elif resolution == "regional":
        df = scale_to_regional_resolution(
            df, region_country_mapping=region_country_mapping, populations=populations
        )
    else:
        raise ValueError("Input resolution is not recognised")
    breakpoint()
    df.T.to_csv(path_to_output, index_label=["id"])
