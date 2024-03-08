import pandas as pd


def create_road_transport_demand_timeseries(
    path_to_input: str, first_year: int, final_year: int, type_name: str,
    power_scaling_factor: float, conversion_factor: float, historic: bool, path_to_output: str
) -> None:
    # Read annual road transport distance into panda dataframe
    df = pd.read_csv(path_to_input, index_col=[0, 1, 2])

    # ASSUME: Every hour of the year the same amount of distance gets driven.
    # Calculate the distance travelled per hour
    series = df["value"] / 8760

    # Process the road transport distance timeseries
    create_timeseries(type_name, power_scaling_factor, conversion_factor, first_year,
                      final_year, series, path_to_output, historic=historic)


def create_timeseries(vehicle: str, power_scaling_factor: float, conversion_factor: float,
                      first_year: int, final_year: int, series: pd.Series,
                      output_path: str, historic: bool):
    ts_index = pd.to_datetime(range(first_year, final_year + 2), format="%Y")
    ts = (
        series
        .xs(vehicle)
        .xs(slice(first_year, final_year + 1), level=1, drop_level=False)
        .unstack("country_code")
        .set_index(ts_index)
        .resample("H")
        .ffill()
        .iloc[:-1] # remove first hour in following year
        .mul(conversion_factor)
        .mul(power_scaling_factor)
        .mul(1 if historic else -1) # historic demand is actually a supply to avoid double counting
    )
    ts_index = pd.to_datetime(ts.index, format="%Y")
    ts = ts.set_index(ts_index)
    ts.to_csv(output_path, index_label="utc-timestamp")


if __name__ == "__main__":
    create_road_transport_demand_timeseries(
        path_to_input=snakemake.input.data,
        power_scaling_factor=snakemake.params.power_scaling_factor,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        type_name=snakemake.params.type_name,
        conversion_factor=snakemake.params.conversion_factor,
        historic=snakemake.params.historic,
        path_to_output=snakemake.output[0],
    )
