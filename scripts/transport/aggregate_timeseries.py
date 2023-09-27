import pandas as pd


def aggregate_timeseries(resolution, timeseries_paths, output_path):
    if resolution == "continental":
        aggregate_continental_timeseries(timeseries_paths, output_path)
    elif resolution == "national":
        aggregate_national_timeseries(timeseries_paths, output_path)
    elif resolution == "regional":
        aggregate_regional_timeseries(timeseries_paths, output_path)


def aggregate_continental_timeseries(timeseries_paths, output_path):
    aggregated_timeseries_df = pd.DataFrame()

    # Loop over all the provided timeseries and read each file into a DataFrame
    for timeseries in timeseries_paths:
        df = pd.read_csv(timeseries, index_col='utc-timestamp', parse_dates=True)

        # Sum the values across all countries for each utc-timestamp
        df['EUR'] = df.sum(axis=1)

        # Extract only the summed 'EUR' column and add it to the aggregated_timeseries_df
        if aggregated_timeseries_df.empty:
            aggregated_timeseries_df = df[['EUR']]
        else:
            aggregated_timeseries_df['EUR'] += df['EUR']

    # Write the created aggregated timeseries DataFrame into a csv file at the desired output path.
    aggregated_timeseries_df.to_csv(output_path, index_label="utc-timestamp")


def aggregate_national_timeseries(timeseries_paths, output_path):
    aggregated_timeseries_df = pd.DataFrame()

    # Loop over all the provided timeseries and read each file into a DataFrame
    for timeseries in timeseries_paths:

        df = pd.read_csv(timeseries, index_col='utc-timestamp', parse_dates=True)

        # If the aggregated timeseries DataFrame is empty, assign the current DataFrame to it
        # Otherwise, add the current DataFrame to the aggregated timeseries DataFrame
        if aggregated_timeseries_df.empty:
            aggregated_timeseries_df = df
        else:
            aggregated_timeseries_df += df
    # Write the created aggregated timeseries dataframe into a csv file at the desired output path.
    aggregated_timeseries_df.to_csv(output_path, index_label="utc-timestamp")


def aggregate_regional_timeseries(timeseries_paths, output_path):
    raise NotImplementedError("Regional road transport (bau) has not yet been implemented!")


if __name__ == "__main__":
    aggregate_timeseries(
        resolution=snakemake.wildcards.resolution,
        timeseries_paths=snakemake.input.electrified_road_transport_timeseries,
        output_path=snakemake.output.road_transport_timeseries)
    aggregate_timeseries(
        resolution=snakemake.wildcards.resolution,
        timeseries_paths=snakemake.input.electrified_road_bau_transport_timeseries,
        output_path=snakemake.output.road_transport_bau_timeseries
    )
