"""Preprocess national electricity load time series."""
from datetime import datetime, timezone

import pandas as pd
import numpy as np
import pycountry


def national_load(path_to_raw_load, number_rows_valid, year, path_to_output):
    """Extracts national load time series for all countries in a specified year."""
    load = read_load_profiles(
        path_to_raw_load=path_to_raw_load,
        number_rows_valid=number_rows_valid
    )
    load = filter_national(load)
    load = select_year_and_fill_gaps(load, year)
    load = handle_outliers(load)
    load.to_csv(
        path_to_output,
        header=True
    )


def read_load_profiles(path_to_raw_load, number_rows_valid):
    """Reads national load data and handles outliers."""
    data = pd.read_csv(path_to_raw_load, nrows=number_rows_valid, parse_dates=[3])
    data = data[(data["variable"] == "load")]
    data = remove_entsoe_power_statistic_data_where_possible(data)
    data.drop(["variable", "attribute"], axis=1, inplace=True)
    return data.pivot(columns="region", index="utc_timestamp", values="data")


def select_year_and_fill_gaps(load_df, year):
    """Selects relevant year then fills in all NaNs with data from other years"""

    year = str(year)  # not sure if it comes in as a string or an int, so we make sure of its type here
    missing_data_regions = set(load_df.columns).difference(
        load_df.loc[year].dropna(axis=1).columns.unique()
    )

    for region in missing_data_regions:
        this_region_df = load_df[region].copy()
        this_region_df[this_region_df < 0] = np.nan  # negative values not allowed

        all_missing_timesteps = this_region_df[year][np.isnan(this_region_df[year])].index
        if len(all_missing_timesteps) <= 48:
            continue
        fill_years = []
        # keep going until almost all gaps are filled and we have a frankenstein's monster of a timeseries
        # We're OK with 2 days of gap (48hrs)
        while this_region_df.loc[all_missing_timesteps].isnull().sum() > 48:
            # Plan to take the next year's data
            avail_years = this_region_df.loc[slice(str(int(year) + 1), None)].dropna().index
            # If needed, take the previous year's data
            if avail_years.empty is True:
                avail_years = this_region_df.loc[slice(None, str(int(year) - 1))].dropna().index
            # There might be nothing!
            if avail_years.empty is True:
                print('no available data for {} national demand'.format(region))
                break

            next_avail_year = avail_years.year.unique()[0]
            missing_timesteps = this_region_df[year][np.isnan(this_region_df[year])].index

            # Ignore February 29th if next available year doesn't have that data available
            if pd.to_datetime(year + '-02-29').date() in missing_timesteps.date:
                if not pd.Period(freq='Y', year=next_avail_year).is_leap_year:
                    missing_timesteps = missing_timesteps[
                        missing_timesteps.date != pd.to_datetime(year + '-02-29').date()
                    ]

            new_data = this_region_df.loc[missing_timesteps.map(lambda dt: dt.replace(year=next_avail_year))]
            this_region_df.loc[str(next_avail_year)] = np.nan
            new_data.index = missing_timesteps
            this_region_df.update(new_data)
            fill_years += [str(next_avail_year)]
        else:
            load_df.update(this_region_df.loc[all_missing_timesteps].to_frame(region))
            print(
                'Country {} has {} missing load values, a working dataset was constructed '
                'from year(s) {}. {} missing timesteps will be filled from nearby values'
                .format(region, len(all_missing_timesteps), ','.join(fill_years), 
                        this_region_df.loc[all_missing_timesteps].notnull().sum())
            )

    return load_df[year]


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
