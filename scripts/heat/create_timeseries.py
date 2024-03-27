import pandas as pd


def create_heat_demand_timeseries(
    path_to_annual_demand: str,
    power_scaling_factor: float,
    first_year: int,
    final_year: int,
    historic: bool,
    path_to_output: str,
) -> None:
    # Read annual heat demand into panda dataframe before inputing timeseries profile
    breakpoint()
    annual_demand = (  # TODO do not just sum over all types of end_use heat
        pd.read_csv(path_to_annual_demand, index_col=[0, 1, 2], parse_dates=[0])
        .xs(slice(str(first_year), str(final_year + 1)), level="year", drop_level=False)
        .groupby(level="year")
        .sum()
        .tz_localize("UTC")
        .resample("H")
        .ffill()
        .iloc[:-1]
    )

    # ASSUME flat profiles for heat demand for now
    # ASSUME that heat to electricity demand is 1:1 ratio for now
    df_timeseries = (
        annual_demand.groupby(by=lambda idx: idx.year)
        .transform(lambda x: x / len(x.index))
        .mul(
            1e6 * power_scaling_factor
        )  # convert from TWh to MWh and then to scaling factor
        .mul(
            1 if historic else -1
        )  # historic demand is actually a supply to avoid double counting
        .tz_localize(None)
        .rename_axis("utc-timestamp")
    )

    assert not df_timeseries.isna().any(
        axis=None
    ), "There are NaN values in the timeseries dataframe"

    df_timeseries.to_csv(path_to_output)


if __name__ == "__main__":
    create_heat_demand_timeseries(
        path_to_annual_demand=snakemake.input.annual_demand,
        power_scaling_factor=snakemake.params.power_scaling_factor,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        historic=snakemake.params.historic,
        path_to_output=snakemake.output[0],
    )
