"""
Preprocess national electricity load time series.

Data is loaded and (optionally) gap filled as follows:

- select only data for countries of interest
- remove outlier data
- select only data for year of interest
- interpolate empty data when the data gap is small
- fill larger data gaps with data from nearby years
- If 29th of February is empty and couldn't previously be filled, use data from the 28th of February
- select most complete dataset for each country from list of data sources, based on a priority order

This process will fail if the 'most complete dataset' for any country has any remaining gaps.

"""

import calendar

import pandas as pd
import pycountry
import logging

logger = logging.getLogger(_name_)
logging.basicConfig(filename='snakemake.log', encoding='utf-8', level=logging.DEBUG)


def national_load(
    path_to_raw_load,
    first_year,
    final_year,
    data_quality_config,
    path_to_output,
    countries,
):
    """Extracts national load time series for all countries in a specified year."""

    load = clean_load_data(
        path_to_raw_load, first_year, final_year, data_quality_config, countries
    )
    load.to_csv(path_to_output, header=True)


def clean_load_data(
    path_to_raw_load, first_year, final_year, data_quality_config, countries
):
    data_sources = data_quality_config["data-source-priority-order"]
    raw_load = read_load_profiles(path_to_raw_load, data_sources)
    filtered_load = filter_countries(raw_load, countries)
    filtered_load = filter_outliers(filtered_load, data_quality_config)
    gap_filled_load = pd.concat(
        fill_gaps_per_source(filtered_load, year, data_quality_config, source)
        for source in data_sources
        for year in range(first_year, final_year + 1)
    ).sort_index()
    return get_source_choice_per_country(
        filtered_load[
            (
                filtered_load.index.get_level_values("utc_timestamp")
                >= f"{first_year}-01-01 00:00"
            )
            & (
                filtered_load.index.get_level_values("utc_timestamp")
                <= f"{final_year}-12-31 23:00"
            )
        ],
        gap_filled_load,
        data_sources,
    )


def read_load_profiles(path_to_raw_load, entsoe_priority):
    """Reads national load data."""
    data = pd.read_csv(path_to_raw_load, parse_dates=["utc_timestamp"])
    load_by_attribute = (
        data[(data.variable == "load") & (data.attribute.isin(entsoe_priority))]
        .set_index(["utc_timestamp", "attribute", "region"])
        .loc[:, "data"]
        .unstack("region")
    )
    return load_by_attribute


def filter_countries(load, countries):
    """Filters data to configured workflow spatial scope and renames country codes to ISO3."""
    # UKM == United Kingdom of Great Britain and Northern Ireland.
    # Other subsets of the UK exist in the data, with UKM being the only one to cover the whole country.
    load.rename(columns={"GB_UKM": "GB"}, inplace=True)
    country_codes = {
        pycountry.countries.lookup(country).alpha_2: pycountry.countries.lookup(
            country
        ).alpha_3
        for country in countries
    }

    national = (
        load.loc[:, country_codes.keys()]
        .rename(columns=country_codes)
        .rename_axis(columns="country_code")
    )

    return national


def filter_outliers(load, data_quality_config):
    normed_load = load / load.mean()
    cleaned_load = load.where(
        (
            normed_load
            >= data_quality_config["outlier-data-thresholds"]["relative-to-mean-min"]
        )
        & (
            normed_load
            <= data_quality_config["outlier-data-thresholds"]["relative-to-mean-max"]
        )
    )
    return cleaned_load


def fill_gaps_per_source(all_load, model_year, data_quality_config, source):
    """For each valid data source, fills in all NaNs by interpolation or with data from other years"""
    source_specific_load = all_load.xs(source, level="attribute")
    source_specific_load = _interpolate_gaps(
        source_specific_load, data_quality_config["max-interpolate-timesteps"]
    )
    # If there is some data in given year, fill in the rest from other years
    if model_year in source_specific_load.index.year.unique():
        source_specific_model_year_load = source_specific_load.loc[str(model_year)]
        source_specific_model_year_load = _fill_gaps_from_other_years(
            source_specific_model_year_load,
            source_specific_load,
            data_quality_config,
            source,
            model_year,
        )
    # If there is no data in given year, return a series of NaNs spanning the whole year
    else:
        source_specific_model_year_load = source_specific_load.reindex(
            pd.date_range(
                f"{model_year}-01-01", f"{model_year}-12-31", freq="H", tz="UTC"
            )
        )

    return source_specific_model_year_load.assign(attribute=source).set_index(
        "attribute", append=True
    )


def _fill_gaps_from_other_years(
    model_year_load, all_load, data_quality_config, source, model_year
):
    missing_data_countries = _countries_with_missing_data_in_model_year(model_year_load)

    for country in missing_data_countries:
        updated_country_series, N_missing_timesteps, fill_years = (
            _fill_missing_data_in_country(
                all_load.loc[:, country],
                model_year,
                data_quality_config["acceptable-year-diff-for-gap-filling"],
            )
        )
        if data_quality_config["fill-29th-feb-from-28th"]:
            updated_country_series = _fill_29th_feb(updated_country_series, model_year)

        model_year_load.loc[:, country] = updated_country_series
        logger.info(
            f"Using data source `{source}`, {country} has {N_missing_timesteps} missing load value(s) in {model_year}. "
            f"A working dataset was constructed from year(s) {', '.join(fill_years)} "
            f"with {updated_country_series.isnull().sum()} remaining empty data points."
        )

    return model_year_load


def _fill_missing_data_in_country(
    country_series, model_year, acceptable_year_diff_for_gap_filling
):
    all_missing_timesteps = _get_index_of_missing_data(country_series, model_year)

    fill_years = []
    years_with_data = (
        country_series.dropna().index.year.drop(model_year, errors="ignore").unique()
    )
    allowed_years_with_data = years_with_data[
        (years_with_data >= model_year - acceptable_year_diff_for_gap_filling)
        & (years_with_data <= model_year + acceptable_year_diff_for_gap_filling)
    ]
    order = (allowed_years_with_data.values - model_year).astype(float)
    order[order < 0] -= 0.1  # negative years = older

    preferred_order_of_years_with_data = allowed_years_with_data.sort_values(
        key=lambda x: abs(order)
    )

    _missing_timesteps = all_missing_timesteps.copy()
    for next_avail_year in preferred_order_of_years_with_data:
        if len(_missing_timesteps) == 0:
            break
        _missing_timesteps = _get_index_of_missing_data(country_series, model_year)
        _missing_timesteps = _ignore_feb_29th(
            model_year, next_avail_year, _missing_timesteps
        )
        new_data = country_series.reindex(
            _missing_timesteps.map(
                lambda dt, next_avail_year=next_avail_year: dt.replace(
                    year=next_avail_year
                )
            )
        )
        new_data.index = _missing_timesteps
        country_series.update(new_data)
        fill_years.append(str(next_avail_year))

    return country_series.loc[str(model_year)], len(all_missing_timesteps), fill_years


def _get_index_of_missing_data(series, model_year):
    attempted_interpolation = series.where(series > 0)[str(model_year)]
    return attempted_interpolation[attempted_interpolation.isnull()].index


def _ignore_feb_29th(this_year, next_avail_year, missing_timesteps):
    """Ignore February 29th if next available year doesn't have that data available."""
    if (
        calendar.isleap(this_year)
        and not calendar.isleap(next_avail_year)
        and pd.to_datetime(f"{this_year}-02-29").date() in missing_timesteps.date
    ):
        return missing_timesteps[
            missing_timesteps.date != pd.to_datetime(f"{this_year}-02-29").date()
        ]
    else:
        return missing_timesteps


def _fill_29th_feb(load, year):
    if f"{year}-02-29" in load.index.strftime("%Y-%m-%d"):
        filler = load.loc[f"{year}-02-28"]
        filler.index = filler.index + pd.DateOffset(1)
        return load.fillna(filler)
    else:
        return load


def _interpolate_gaps(load, interpolate_timesteps):
    if interpolate_timesteps == 0:
        return load
    else:
        return load.interpolate(limit=interpolate_timesteps, limit_direction="both")


def _countries_with_missing_data_in_model_year(data):
    model_year_missing_data = data.where(data > 0).isnull().sum() > 0
    return model_year_missing_data[model_year_missing_data].index


def get_source_choice_per_country(raw_load, gap_filled_load, entsoe_priority):
    """
    OPSD collects load data from different data sources.
    The hourly data can vary between these countries, with some proving more reliable than others.
    See https://nbviewer.jupyter.org/github/Open-Power-System-Data/datapackage_timeseries/blob/2020-10-06/main.ipynb
    """

    source_choice = (
        gap_filled_load.notnull()
        .groupby(level="attribute")
        .sum()
        .loc[entsoe_priority]
        .idxmax()
        .to_frame("attribute")
        .set_index("attribute", append=True)
    )

    print(
        "Using the following data sources for national load:\n{}".format(
            "\n".join([f"{idx[0]}: {idx[1]}" for idx in source_choice.index])
        )
    )
    new_load = _select_load_by_source_priority(gap_filled_load, source_choice)
    raw_load_for_comparison = _select_load_by_source_priority(raw_load, source_choice)

    if new_load.isnull().any().any():
        bad_index_values = new_load.isnull().stack()
        error_msg = "Gap filling thresholds do not allow for a complete load dataset to be produced. "
        if bad_index_values[bad_index_values].size < 100:
            error_msg = (
                error_msg
                + f"Remaining empty data: {bad_index_values[bad_index_values].index.to_list()}"
            )
        raise AssertionError(error_msg)

    else:
        print(
            "Gap filling methods lead to the following relative increase in output load compared to input data\n"
            f"{new_load.sum() / raw_load_for_comparison.sum()}"
        )
    return new_load


def _select_load_by_source_priority(load, source_priority):
    load_filtered_by_priority_source = (
        load.unstack("attribute")
        .loc[:, source_priority.index]
        .droplevel("attribute", axis="columns")
    )
    assert load_filtered_by_priority_source.columns.duplicated().sum() == 0

    return load_filtered_by_priority_source


if __name__ == "__main__":
    national_load(
        path_to_raw_load=snakemake.input.load,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        data_quality_config=snakemake.params.data_quality_config,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0],
    )
