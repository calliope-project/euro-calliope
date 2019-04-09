"""Preprocess national electricity load time series."""
from datetime import datetime, timezone

import pandas as pd
import numpy as np
import pycountry


def national_load(path_to_raw_load, number_rows_valid, year, path_to_output):
    """Extracts national load time series for all countries in a specified year."""
    load = read_load_profiles(
        path_to_raw_load=path_to_raw_load,
        number_rows_valid=number_rows_valid,
        start=datetime(year, 1, 1, tzinfo=timezone.utc),
        end=datetime(year + 1, 1, 1, tzinfo=timezone.utc)
    )
    load = filter_national(load)
    check_completeness(load)
    load = handle_outliers(load)
    load.to_csv(
        path_to_output,
        header=True
    )


def read_load_profiles(path_to_raw_load, number_rows_valid, start, end):
    """Reads national load data and handles outliers."""
    data = pd.read_csv(path_to_raw_load, nrows=number_rows_valid, parse_dates=[3])
    data = data[(data["variable"] == "load")]
    data = data[(data.utc_timestamp >= start) &
                (data.utc_timestamp < end)]
    data = remove_entsoe_power_statistic_data_where_possible(data)
    data.drop(["variable", "attribute"], axis=1, inplace=True)
    return data.pivot(columns="region", index="utc_timestamp", values="data")


def remove_entsoe_power_statistic_data_where_possible(load):
    sorted_load = load.sort_values(
        "attribute",
        ascending=False
    ) # will end with entsoe-transparency ahead of entsoe-power-statistics
    return sorted_load.drop_duplicates(["region", "utc_timestamp"], keep="first")


def filter_national(load):
    countries = [iso2 for iso2 in load.columns.unique() if len(iso2) == 2]
    national = load.loc[:, countries].copy()
    national.columns.name = "country_code"
    return national.rename(columns=lambda iso2: pycountry.countries.lookup(iso2).alpha_3)


def check_completeness(load):
    for country in load.columns[load.isnull().any()].tolist():
        print("Country {} has missing load values.".format(pycountry.countries.lookup(country).name))


def handle_outliers(all_time_series):
    # considers all data < 0.25 * mean and > 2 * mean invalid and replaces with last valid value
    normed_load = all_time_series / all_time_series.mean()
    all_time_series[(normed_load < 0.25) | (normed_load > 2)] = np.nan
    return all_time_series.fillna(method="ffill")


if __name__ == "__main__":
    national_load(
        path_to_raw_load=snakemake.input.load[0],
        number_rows_valid=snakemake.params.number_rows_valid,
        year=snakemake.params.year,
        path_to_output=snakemake.output[0]
    )
