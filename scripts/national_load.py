"""
Preprocess national electricity load time series.

Data is loaded and (optionally) gap filled as follows:

- select only data for countries of interest
- remove outlier data
- select only data for year of interest
- interpolate empty data when the data gap is small
- fill larger data gaps with data from nearby years
- If 29th of February is empty and couldn't previously be filled, use data from the 28th of February
- interpolate across remaining data gaps if they are now small enough
- select most complete dataset for each country from list of data sources, based on a priority order

This process will fail if the 'most complete dataset' for any country has any remaining gaps.

"""
import calendar

import pandas as pd
import pycountry


def national_load(
    path_to_raw_load, year, data_quality_config, path_to_output, countries
):
    """Extracts national load time series for all countries in a specified year."""
    raw_load = read_load_profiles(path_to_raw_load, data_quality_config["data-source-priority-order"])
    filtered_load = filter_countries(raw_load, countries)
    filtered_load = filter_outliers(filtered_load, data_quality_config)
    gap_filled_load = fill_gaps_per_source(filtered_load, year, data_quality_config)
    load = get_source_choice_per_country(
        filtered_load.loc[str(year)],
        gap_filled_load,
        data_quality_config["data-source-priority-order"]
    )

    load.to_csv(path_to_output, header=True)


def read_load_profiles(path_to_raw_load, entsoe_priority):
    """Reads national load data and handles outliers."""
    data = pd.read_csv(path_to_raw_load, parse_dates=["utc_timestamp"])
    load_by_attribute = (
        data
        [(data.variable == "load") & (data.attribute.isin(entsoe_priority))]
        .set_index(["utc_timestamp", "attribute", "region"])
        .loc[:, "data"]
        .unstack("region")
    )
    return load_by_attribute


def filter_countries(load, countries):
    # UKM == United Kingdom of Great Britain and Northern Ireland.
    # Other subsets of the UK exist in the data, with this being the only one to cover the whole country.
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


def fill_gaps_per_source(all_load, model_year, data_quality_config):
    """For each valid data source, fills in all NaNs by interpolation or with data from other years"""
    source_loads = []
    for source in data_quality_config["data-source-priority-order"]:
        source_specific_load = all_load.xs(source, level="attribute")
        source_specific_load = _interpolate_gaps(source_specific_load, data_quality_config["interpolate-hours"])
        source_specific_model_year_load = source_specific_load.loc[str(model_year)]
        source_specific_model_year_load = _fill_gaps_from_other_years(
            source_specific_model_year_load, source_specific_load, data_quality_config, source, model_year
        )
        if data_quality_config["fill-29th-feb-from-28th"]:
            source_specific_model_year_load = _fill_29th_feb(source_specific_model_year_load, model_year)

        source_loads.append(
            _interpolate_gaps(source_specific_model_year_load, data_quality_config["interpolate-hours"])
            .assign(attribute=source)
            .set_index("attribute", append=True)
        )

    return pd.concat(source_loads)


def _fill_gaps_from_other_years(model_year_load, all_load, data_quality_config, source, model_year):
    missing_data_countries = _countries_with_missing_data_in_model_year(model_year_load)

    for country in missing_data_countries:
        updated_country_series, N_missing_timesteps, fill_years = _fill_missing_data_in_country(
            all_load.loc[:, country], model_year, data_quality_config["acceptable-year-diff-for-gap-filling"]
        )
        model_year_load.loc[:, country] = updated_country_series
        print(
            f"Using data source `{source}`, {country} has {N_missing_timesteps} missing load value(s). "
            f"A working dataset was constructed from year(s) {', '.join(fill_years)} "
            f"with {updated_country_series.isnull().sum()} remaining empty data points."
        )
    return model_year_load


def _fill_missing_data_in_country(country_series, model_year, acceptable_year_diff_for_gap_filling):
    all_missing_timesteps = _get_index_of_missing_data(country_series, model_year)

    fill_years = []
    years_with_data = country_series.dropna().index.year.drop(model_year, errors="ignore").unique()
    allowed_years_with_data = years_with_data[
        (years_with_data >= model_year - acceptable_year_diff_for_gap_filling) &
        (years_with_data <= model_year + acceptable_year_diff_for_gap_filling)
    ]
    def __give_older_years_lower_priority(years):
        order = (years.values - model_year).astype(float)
        order[order < 0] -= 0.1 # negative years = older
        return abs(order)

    preferred_order_of_years_with_data = list(allowed_years_with_data.sort_values(key=__give_older_years_lower_priority))

    _missing_timesteps = all_missing_timesteps.copy()
    for next_avail_year in preferred_order_of_years_with_data:
        if len(_missing_timesteps) == 0:
            break
        _missing_timesteps = _get_index_of_missing_data(country_series, model_year)
        _missing_timesteps = _ignore_feb_29th(model_year, next_avail_year, _missing_timesteps)
        new_data = country_series.reindex(_missing_timesteps.map(lambda dt: dt.replace(year=next_avail_year)))
        new_data.index = _missing_timesteps
        country_series.update(new_data)
        fill_years += [str(next_avail_year)]

    return country_series.loc[str(model_year)], len(all_missing_timesteps), fill_years


def _get_index_of_missing_data(series, model_year):
    attempted_interpolation = series.where(series >= 0)[str(model_year)]
    return attempted_interpolation[attempted_interpolation.isnull()].index


def _ignore_feb_29th(this_year, next_avail_year, missing_timesteps):
    """ Ignore February 29th if next available year doesn't have that data available."""
    if (calendar.isleap(this_year)
        and not calendar.isleap(next_avail_year)
        and pd.to_datetime(f'{this_year}-02-29').date() in missing_timesteps.date
    ):
        return missing_timesteps[
            missing_timesteps.date != pd.to_datetime(f'{this_year}-02-29').date()
        ]
    else:
        return missing_timesteps


def _fill_29th_feb(load, year):
    if pd.to_datetime(f'{year}-02-29').date() in load.index.date:
        filler = load.loc[f'{year}-02-28']
        to_fill = load.loc[f'{year}-02-29']
        filler.index = to_fill.index
        load.loc[f'{year}-02-29'] = to_fill.fillna(filler)
    return load


def _interpolate_gaps(load, interpolate_hours):
    return load.interpolate(limit=interpolate_hours)


def _countries_with_missing_data_in_model_year(data):
    model_year_missing_data = data.where(data >= 0).isnull().sum() > 0
    return model_year_missing_data[model_year_missing_data].index


def filter_outliers(load, data_quality_config):
    normed_load = load / load.mean()
    cleaned_load = load.where(
        (normed_load >= data_quality_config["outlier-data-thresholds"]["relative-to-mean-min"]) &
        (normed_load <= data_quality_config["outlier-data-thresholds"]["relative-to-mean-max"])
    )
    return cleaned_load


def get_source_choice_per_country(raw_load, gap_filled_load, entsoe_priority):
    """
    OPSD collects load data from different data sources.
    The hourly data can vary between these countries, with some proving more reliable than others.
    See https://nbviewer.jupyter.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb
    """

    source_choice = (
        gap_filled_load
        .notnull()
        .sum(level="attribute")
        .loc[entsoe_priority]
        .idxmax()
        .to_frame("attribute")
        .set_index("attribute", append=True)
    )

    print(
        "Using the following data sources for national load:\n{}"
        .format("\n".join([f"{idx[0]}: {idx[1]}" for idx in source_choice.index]))
    )
    new_load = _select_load_by_source_priority(gap_filled_load, source_choice)
    raw_load_for_comparison = _select_load_by_source_priority(raw_load, source_choice)

    if new_load.isnull().sum().all() == 0:
        print(
            "Gap filling methods lead to the following relative differences between input and output data\n"
            f"{new_load.sum() / raw_load_for_comparison.sum()}"
        )
    else:
        raise AssertionError(
            "Gap filling thresholds do not allow for a complete load dataset to be produced. "
            f"Remaining empty data: {new_load[new_load.isnull()]}"
        )
    return new_load


def _select_load_by_source_priority(load, source_priority):
    return (
        load
        .stack()
        .align(source_priority)[0]
        .unstack('country_code')
        .droplevel("attribute")
    )


if __name__ == "__main__":
    national_load(
        path_to_raw_load=snakemake.input.load,
        year=snakemake.params.year,
        data_quality_config=snakemake.params.data_quality_config,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0]
    )
