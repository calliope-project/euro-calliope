import pandas as pd


def aggregate_timeseries(resolution: str, paths_to_input: list[str], path_to_output: str) -> None:
    if resolution == "continental":
        ts = aggregate_continental_timeseries(paths_to_input)
    elif resolution == "national":
        ts = aggregate_national_timeseries(paths_to_input)
    elif resolution == "regional":
        ts = aggregate_regional_timeseries(paths_to_input)

    ts.to_csv(path_to_output)


def aggregate_continental_timeseries(paths_to_input: list[str]) -> pd.DataFrame:
    ts = aggregate_national_timeseries(paths_to_input)
    return ts.sum(axis=1).rename("EUR")


def aggregate_national_timeseries(paths_to_input: list[str]) -> pd.DataFrame:
    all_ts = [
        pd.read_csv(path, index_col='utc-timestamp', parse_dates=True)
        for path in paths_to_input
    ]
    return sum(all_ts)


def aggregate_regional_timeseries(paths_to_input: list[str]) -> pd.DataFrame:
    # TODO implement
    raise NotImplementedError("Regional road transport (historic) has not yet been implemented!")


if __name__ == "__main__":
    aggregate_timeseries(
        resolution=snakemake.wildcards.resolution,
        paths_to_input=snakemake.input.time_series,
        path_to_output=snakemake.output[0]
    )
