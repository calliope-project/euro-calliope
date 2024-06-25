import pandas as pd
import pycountry


def create_road_transport_demand_timeseries(
    path_to_annual_data: str,
    path_to_timeseries: str,
    first_year: int,
    final_year: int,
    vehicle_type: str,
    country_neighbour_dict: dict[str, list[str]],
    power_scaling_factor: float,
    conversion_factor: float,
    country_codes: list[str],
    path_to_output: str,
) -> None:
    """This function reads historical electrified transport annual data and creates a timeseries to act as a supply"""

    # Read annual road transport distance into panda dataframe

    df_annual = (
        pd.read_csv(path_to_annual_data, index_col=[0, 1, 2], parse_dates=[2])
        .squeeze()
        .xs(vehicle_type)
        .xs(slice(str(first_year), str(final_year + 1)), level="year", drop_level=False)
        .unstack("country_code")
        .tz_localize("UTC")
        .resample("H")
        .ffill()
        .iloc[:-1]
    )

    if vehicle_type in ["passenger-cars", "motorcycles", "light-duty-vehicles"]:
        # Use RAMP time series profiles for small and light vehicles.
        df_timeseries = (
            pd.read_csv(path_to_timeseries, index_col=[0, 1, 2], parse_dates=[0])
            .xs(slice(first_year, final_year), level="year")
            .unstack("country_code")
            .droplevel(level=0, axis="columns")
            .groupby(by=lambda idx: idx.year)
            .transform(lambda x: x / x.sum())
            .pipe(fill_empty_country, country_neighbour_dict)
            .tz_localize("UTC")
            .mul(df_annual)
        )
    elif vehicle_type in ["heavy-duty-vehicles", "coaches-and-buses"]:
        # ASSUME flat profiles for heavy transport.
        df_timeseries = df_annual.groupby(by=lambda idx: idx.year).transform(
            lambda x: x / len(x.index)
        )

    else:
        raise ValueError(f"vehicle_type {vehicle_type} is not supported")

    df_timeseries = (
        df_timeseries.mul(conversion_factor)
        .mul(power_scaling_factor)
        .loc[:, country_codes]
        .tz_localize(None)
        .rename_axis("utc-timestamp")
    )

    assert not df_timeseries.isna().any(
        axis=None
    ), "There are NaN values in the timeseries dataframe"

    df_timeseries.to_csv(path_to_output)


def fill_empty_country(df, country_neighbour_dict):
    for country, neighbours in country_neighbour_dict.items():
        if country in df.columns:
            print(f"Country {country} is already in dataframe")
            continue
        print(f"Country {country} is not in dataframe, filling with neighbours")
        df[country] = df[neighbours].mean(axis=1)
    return df


if __name__ == "__main__":
    create_road_transport_demand_timeseries(
        path_to_annual_data=snakemake.input.annual_data,
        path_to_timeseries=snakemake.input.timeseries,
        power_scaling_factor=snakemake.params.power_scaling_factor,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        vehicle_type=snakemake.wildcards.vehicle_type,
        conversion_factor=snakemake.params.conversion_factor,
        path_to_output=snakemake.output[0],
        country_codes=[
            pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
        ],
        country_neighbour_dict=snakemake.params.country_neighbour_dict,
    )
