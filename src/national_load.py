"""Preprocess national electricity load time series."""
from datetime import datetime, timezone

import pandas as pd
import numpy as np
import pycountry

SOURCE_PRIORITY = [
    'actual_entsoe_power_statistics',
    'actual_entsoe_transparency',
    'actual_tso',
    'actual_net_consumption_tso'
]

def national_load(path_to_raw_load, number_rows_valid, year, path_to_output):
    """Extracts national load time series for all countries in a specified year."""
    load = read_load_profiles(
        path_to_raw_load=path_to_raw_load,
        number_rows_valid=number_rows_valid,
        start=datetime(year, 1, 1, tzinfo=timezone.utc),
        end=datetime(year + 1, 1, 1, tzinfo=timezone.utc)
    )
    if year < 2017: # data for Albania before 2017 is missing
        load["AL"] = read_albania(path_to_raw_load, number_rows_valid, other_ts=load)
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
    data = select_statistics_by_source_priority(data)
    return data.unstack("region")


def read_albania(path_to_raw_load, number_rows_valid, other_ts):
    albania = read_load_profiles(
        path_to_raw_load="data/automatic/raw-load-data.csv",
        number_rows_valid=20950674,
        start=datetime(2017, 1, 1, tzinfo=timezone.utc),
        end=datetime(2017 + 1, 1, 1, tzinfo=timezone.utc)
    )["AL"]
    if len(other_ts.index) == 8784: # leap year
        albania = pd.concat([albania[:"2017-02"], albania["2017-02-28"], albania["2017-03":]])
    albania.index = other_ts.index
    return albania


def select_statistics_by_source_priority(load):
    """
    Choosing `entsoe_power_statistics` as main source since OPSD states:
        The two sources differ Values on PS (~500 TWh annaually in Germany) are
        usually slightly higher than on the TP (~490 TWh). The reason probably
        lies with different reporting deadlines: Values on the TP have to be
        reported "no later than one hour after the end of the operating period".
        For the PS, the data is published with a delay of up to 3 months,
        which might allow for more accurate metering.
        For a comparison of the two sources see Hirth, et al. (2018).
    See https://nbviewer.jupyter.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb for more info.
    """
    load_by_attribute = (
        load
        .set_index(["region", "utc_timestamp", "attribute"])
        ["data"]
        .unstack("attribute")
    )
    load_top_priority = load_by_attribute[SOURCE_PRIORITY[0]]
    for source in SOURCE_PRIORITY[1:]:
        load_top_priority = load_top_priority.fillna(load_by_attribute[source])

    return load_top_priority


def filter_national(load):
    load.rename(columns={"GB_UKM": "GB"}, inplace=True)
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
