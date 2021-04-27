"""Preprocess national electricity load time series."""
import calendar

import pandas as pd
import numpy as np
import pycountry


def national_load(
    path_to_raw_load, entsoe_priority, year, acceptable_gap_hours,
    outlier_thresholds, path_to_output, countries
):
    """Extracts national load time series for all countries in a specified year."""
    load = read_load_profiles(path_to_raw_load, entsoe_priority)
    load = filter_national(load, countries)
    load = select_year_and_fill_gaps(load, year, acceptable_gap_hours)
    load = handle_outliers(load, outlier_thresholds)
    load.to_csv(path_to_output, header=True)


def read_load_profiles(path_to_raw_load, entsoe_priority):
    """Reads national load data and handles outliers."""
    data = pd.read_csv(path_to_raw_load, parse_dates=["utc_timestamp"])
    data = data[(data["variable"] == "load")]
    data = select_statistics_by_source_priority(data, entsoe_priority)
    return data.unstack("region")


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
        updated_country_series, N_missing_timesteps, fill_years = fill_missing_data_in_country(
            load_df.loc[:, country], model_year, acceptable_gap_hours
        )
        load_df.update(updated_country_series.to_frame(country))
        print(
            f"Country {country} has {N_missing_timesteps} missing load value(s). "
            f"A working dataset was constructed from year(s) {', '.join(fill_years)} "
            f"with {updated_country_series.isnull().sum()} remaining empty data points."
        )

    return load_df.loc[str(model_year)]


def _get_index_of_missing_data(series, model_year):
    return series[(series < 0) | series.isna()][str(model_year)].index


def _ignore_feb_29th(this_year, next_avail_year, missing_timesteps):
    """ Ignore February 29th if next available year doesn't have that data available """
    if (calendar.isleap(this_year)
        and not calendar.isleap(next_avail_year)
        and pd.to_datetime(f'{this_year}-02-29').date() in missing_timesteps.date
    ):
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

    def _give_older_years_lower_priority(years):
        year_order = []
        for year in years:
            if year - model_year < 0:
                year_order.append(abs(year - model_year) + 0.1)
            else:
                year_order.append(year - model_year)
        return year_order
    preferred_order_of_years_with_data = list(years_with_data.sort_values(key=_give_older_years_lower_priority, ascending=False))

    _missing_timesteps = all_missing_timesteps.copy()
    while len(_missing_timesteps) > acceptable_gap_hours and preferred_order_of_years_with_data != []:
        next_avail_year = preferred_order_of_years_with_data.pop(0)
        _missing_timesteps = _get_index_of_missing_data(country_series, model_year)
        _missing_timesteps = _ignore_feb_29th(model_year, next_avail_year, _missing_timesteps)
        new_data = country_series.reindex(_missing_timesteps.map(lambda dt: dt.replace(year=next_avail_year)))
        new_data.index = _missing_timesteps
        country_series.update(new_data)
        fill_years += [str(next_avail_year)]

    assert country_series.loc[str(model_year)].isna().sum() <= acceptable_gap_hours, f"Not enough load data for country {country_series.name}"

    return country_series.loc[all_missing_timesteps], len(all_missing_timesteps), fill_years


def select_statistics_by_source_priority(load, entsoe_priority):
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
        .loc[:, "data"]
        .unstack("attribute")
    )
    load_top_priority = load_by_attribute[entsoe_priority[0]]
    if len(entsoe_priority) > 1:
        for source in entsoe_priority[1:]:
            load_top_priority = load_top_priority.fillna(load_by_attribute[source])

    return load_top_priority


def filter_national(load, countries):
    load.rename(columns={"GB_UKM": "GB"}, inplace=True)
    country_codes = {
        pycountry.countries.lookup(country).alpha_2: pycountry.countries.lookup(country).alpha_3
        for country in countries
    }

    national = (
        load
        .loc[:, country_codes.keys()]
        .rename(columns=country_codes)
        .rename_axis(columns="country_code")
    )

    return national


def handle_outliers(all_time_series, outlier_thresholds):
    # considers all data < 0.25 * mean and > 2 * mean invalid and replaces with last valid value
    normed_load = all_time_series / all_time_series.mean()
    all_time_series[(normed_load < outlier_thresholds["relative-to-mean-min"]) | (normed_load > outlier_thresholds["relative-to-mean-max"])] = np.nan
    # check that this outlier handling won't vary any country's annual load by more than a percentage deviation
    assert (
        abs((all_time_series.interpolate().sum() - all_time_series.sum()) / all_time_series.sum()) * 100
        <= outlier_thresholds["percentage-deviation-post-cleaning-max"]
    ).all()
    return all_time_series.interpolate()


if __name__ == "__main__":
    national_load(
        path_to_raw_load=snakemake.input.load,
        entsoe_priority=snakemake.params.entsoe_priority,
        year=snakemake.params.year,
        acceptable_gap_hours=snakemake.params.acceptable_gap_hours,
        outlier_thresholds=snakemake.params.outlier_thresholds,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
