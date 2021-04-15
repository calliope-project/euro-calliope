"""Preprocess national electricity load time series."""
from datetime import datetime, timezone

import pandas as pd
import numpy as np
import pycountry


def national_load(path_to_raw_load, year, acceptable_gap_hours, outlier_thresholds, path_to_output):
    """Extracts national load time series for all countries in a specified year."""
    load = read_load_profiles(path_to_raw_load)
    load = filter_national(load)
    load = select_year_and_fill_gaps(load, year, acceptable_gap_hours)
    load = handle_outliers(load, outlier_thresholds)
    load.to_csv(path_to_output, header=True)


def read_load_profiles(path_to_raw_load):
    """Reads national load data and handles outliers."""
    data = pd.read_csv(path_to_raw_load, parse_dates=["utc_timestamp"])
    data = data[(data["variable"] == "load")]
    data = remove_entsoe_power_statistic_data_where_possible(data)
    data.drop(["variable", "attribute"], axis=1, inplace=True)
    return data.pivot(columns="region", index="utc_timestamp", values="data")


def columns_with_missing_data_in_model_year(data, model_year, acceptable_gap_hours):
    model_year_missing_data = data.where(data >= 0).loc[str(model_year)]
    return model_year_missing_data.dropna(
        axis=1, thresh=len(model_year_missing_data) - acceptable_gap_hours
    )


def select_year_and_fill_gaps(load_df, model_year, acceptable_gap_hours):
    """Selects relevant year then fills in all NaNs with data from other years"""
    model_year_missing_data = columns_with_missing_data_in_model_year(
        load_df, model_year, acceptable_gap_hours
    )
    missing_data_countries = set(load_df.columns).difference(model_year_missing_data.columns)

    for country in missing_data_countries:
        try:
            updated_country_series, N_missing_timesteps, fill_years = fill_missing_data_in_country(
                load_df.loc[:, country], model_year, acceptable_gap_hours
            )
        except AssertionError:
            print(f"Not enough load data for country {country}")
        else:
            load_df.update(updated_country_series.to_frame(country))
            print(
                f"Country {country} has {N_missing_timesteps} missing load value(s). "
                f"A working dataset was constructed from year(s) {', '.join(fill_years)} "
                f"with {updated_country_series.isnull().sum()} remaining empty data points."
            )

    return load_df.loc[str(model_year)]


def _get_index_of_missing_data(series, model_year):
    return series[(series < 0) | series.isna()][str(model_year)].index


def _handle_feb_29th(this_year, next_avail_year, missing_timesteps):
    """ Ignore February 29th if next available year doesn't have that data available """
    if (pd.Period(freq='Y', year=this_year).is_leap_year and
            not pd.Period(freq='Y', year=next_avail_year).is_leap_year and
            pd.to_datetime(f'{this_year}-02-29').date() in missing_timesteps.date):
        return missing_timesteps[
            missing_timesteps.date != pd.to_datetime(f'{this_year}-02-29').date()
        ]
    else:
        return missing_timesteps


def fill_missing_data_in_country(country_series, model_year, acceptable_gap_hours):
    all_missing_timesteps = _get_index_of_missing_data(country_series, model_year)

    fill_years = []
    # keep going until almost all gaps are filled and we have a frankenstein's monster of a timeseries
    # We're OK with some gap, as defined by 'acceptable_gap_hours'
    years_with_data = country_series.dropna().index.year.drop(model_year, errors="ignore").unique()

    def _get_year_order(years):
        year_order = []
        for year in years:
            if year - model_year < 0:
                year_order.append(abs(year - model_year) + 0.1)
            else:
                year_order.append(year - model_year)
        return year_order
    preferred_order_of_years_with_data = list(years_with_data.sort_values(key=_get_year_order, ascending=False))

    while country_series.loc[all_missing_timesteps].isnull().sum() > acceptable_gap_hours:
        if preferred_order_of_years_with_data == []:
            break

        next_avail_year = preferred_order_of_years_with_data.pop(0)
        _missing_timesteps = _get_index_of_missing_data(country_series, model_year)
        _missing_timesteps = _handle_feb_29th(model_year, next_avail_year, _missing_timesteps)
        new_data = country_series.reindex(_missing_timesteps.map(lambda dt: dt.replace(year=next_avail_year)))
        new_data.index = _missing_timesteps
        country_series.update(new_data)
        fill_years += [str(next_avail_year)]

    assert columns_with_missing_data_in_model_year(
        country_series.to_frame(), model_year, acceptable_gap_hours
    ).empty

    return country_series.loc[all_missing_timesteps], len(all_missing_timesteps), fill_years

def remove_entsoe_power_statistic_data_where_possible(load):
    sorted_load = load.sort_values(
        "attribute",
        ascending=False
    ) # will end with entsoe-transparency ahead of entsoe-power-statistics
    return sorted_load.drop_duplicates(["region", "utc_timestamp"], keep="first")


def filter_national(load):
    load.rename(columns={"GB_UKM": "GB"}, inplace=True)
    countries = [iso2 for iso2 in load.columns.unique() if iso2 in [i.alpha_2 for i in pycountry.countries]]
    national = load.loc[:, countries].copy()
    national.columns.name = "country_code"
    return national.rename(columns=lambda iso2: pycountry.countries.lookup(iso2).alpha_3)


def handle_outliers(all_time_series, outlier_thresholds):
    # considers all data < 0.25 * mean and > 2 * mean invalid and replaces with last valid value
    normed_load = all_time_series / all_time_series.mean()
    all_time_series[(normed_load < outlier_thresholds["relative-to-mean-min"]) | (normed_load > outlier_thresholds["relative-to-mean-max"])] = np.nan
    # check that this outlier handling won't vary any country's annual load by more than +/- 1%
    assert (abs((all_time_series.interpolate().sum() - all_time_series.sum()) / all_time_series.sum()) <= outlier_thresholds["abolute-deviation-post-cleaning-max"]).all()
    return all_time_series.interpolate()


if __name__ == "__main__":
    national_load(
        path_to_raw_load=snakemake.input.load[0],
        year=snakemake.params.year,
        acceptable_gap_hours=snakemake.params.acceptable_gap_hours,
        outlier_thresholds=snakemake.params.outlier_thresholds,
        path_to_output=snakemake.output[0]
    )
