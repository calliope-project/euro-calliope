"""Utility functions."""
import pycountry
import pandas as pd
import xarray as xr

from string import digits


def convert_country_code(input_country, output="alpha3"):
    """
    Converts input country code or name into either either a 2- or 3-letter code.

    ISO alpha2: alpha2
    ISO alpha2 with Eurostat codes: alpha2_eurostat
    ISO alpha3: alpha3

    """

    if input_country.lower() == "el":
        input_country = "gr"
    elif input_country.lower() == "uk":
        input_country = "gb"
    elif (
        input_country.lower() == "bh"
    ):  # this is a weird country code used in the biofuels dataset
        input_country = "ba"

    if output == "alpha2":
        return pycountry.countries.lookup(input_country).alpha_2

    if output == "alpha2_eurostat":
        result = pycountry.countries.lookup(input_country).alpha_2
        if result == "GB":
            return "UK"
        elif result == "GR":
            return "EL"
        else:
            return result

    if output == "alpha3":
        return pycountry.countries.lookup(input_country).alpha_3


def convert_valid_countries(country_codes: list, output: str = "alpha3") -> dict:
    """
    Convert a list of country codes / names to a list of uniform ISO coded country
    codes. If an input item isn't a valid country (e.g. "EU27") then print the code and
    continue, instead of raising an exception

    Args:
        country_codes (list):
            Strings defining country codes / names
            (["France", "FRA", "FR"] will all be treated the same)

    Returns:
        dict: Mapping from input country code/name to output country code for all valid input countries
    """

    mapped_codes = {}
    for country_code in country_codes:
        try:
            mapped_codes[country_code] = convert_country_code(country_code, output=output)
        except LookupError:
            print(f"Skipping country/region {country_code} in annual energy balances")
            continue
    return mapped_codes


def rename_and_groupby(
    da: xr.DataArray,
    rename_dict: dict,
    dim_name: str,
    new_dim_name: str = None,
    dropna: bool = False,
    keep_non_renamed: bool = False,
) -> xr.DataArray:
    """
    Take an xarray dataarray and rename the contents of a given dimension
    as well as (optionally) rename that dimension.
    If renaming the contents has some overlap (e.g. {'foo' : 'A', 'bar': 'A'})
    then the returned dataarray will be grouped over the new dimension items
    (by summing the data).

    Args:
        da (xr.DataArray):
            Input dataarray with the dimension "dim_name".
        rename_dict (dict):
            Dictionary to map items in the dimension "dim_name" to new names ({"old_item_name": "new_item_name"}).
        dim_name (str):
            Dimension on which to rename items.
        new_dim_name (str, optional): Defaults to None.
            If not None, rename the dimension "dim_name" to the given string.
        dropna (bool, optional): Defaults to False.
            If True, drop any items in "dim_name" after renaming/grouping which have all NaN values along all other dimensions.
        keep_non_renamed (bool, optional): Defaults to False.
            If False, any item in "dim_name" that is not referred to in "rename_dict" will be removed from that dimension in the returned array.
    Returns:
        (xr.DataArray): Same as "da" but with the items in "dim_name" renamed and possibly a. grouped, b. "dim_name" itself renamed.
    """
    rename_series = pd.Series(rename_dict).rename_axis(index=dim_name)
    if keep_non_renamed is True:
        existing_dim_items = da[dim_name].to_series()
        rename_series = rename_series.reindex(existing_dim_items).fillna(existing_dim_items)

    if new_dim_name is None:
        new_dim_name = f"_{dim_name}"  # placeholder that we'll revert
        revert_dim_name = True
    else:
        revert_dim_name = False

    rename_da = xr.DataArray(rename_series.rename(new_dim_name))
    da = (
        da
        .reindex({dim_name: rename_da[dim_name]})
        .groupby(rename_da)
        .sum(dim_name, skipna=True, min_count=1, keep_attrs=True)
    )
    if revert_dim_name:
        da = da.rename({new_dim_name: dim_name})
        new_dim_name = dim_name
    if dropna:
        da = da.dropna(new_dim_name, how="all")
    return da


def merge_da(da_list: list, merged_da_name: str = None) -> xr.DataArray:
    """
    Merge dataArrays with the same dimensions but different dimension items
    into a single xarray datarray

    Args:
        da_list (list): list of xarray dataArrays
        merged_da_name (str, optional): Defaults to None.
            Name of merged datarray

    Returns:
        xr.DataArray:
            Merged Datarray, in which all dimensions contain all items defined in the
            arrays in `da_list`

    """
    datasets = [da.rename("var") for da in da_list]
    return xr.merge(datasets, combine_attrs="no_conflicts")["var"].rename(merged_da_name)


def to_numeric(series):
    """
    Clean up a pandas.Series which was parsed as strings, but is really numeric:

    1. replace "-" for "NaN" into numbers and NaNs
    2. removes random superscript attached to numbers
       (e.g. pointing to footnotes in an excel), "1000c" -> 1000

    Returns a numeric pandas.Series.

    """
    series = series.astype(str).str.extract("(\\-*\\d+\\.*\\d*)")[0]
    return pd.to_numeric(series, errors="coerce")


def gwh_to_tj(array):
    """Convert GWh to TJ"""
    return array * 3.6


def gwh_to_twh(array):
    """Convert GWh to TWh"""
    return array / 1000


def pj_to_twh(array):
    """Convert PJ to TWh"""
    return array / 3.6


def tj_to_twh(array):
    """Convert TJ to TWh"""
    return pj_to_twh(array) / 1000


def ktoe_to_twh(array):
    """Convert KTOE to TWH"""
    return array * 1.163e-2


def read_eurostat_tsv(path_to_tsv, slice_idx=None, slice_lvl=None):
    """

    Read a typical tab-delimited file from EUROSTAT. These have a specific structure
    where the data is tab-delimited but the multi-index data is comma delimited.
    This function also prepares the data in the expectation that it is all numeric
    and that the columns are given as years (a standard EUROSTAT format)

    Parameters
    ---------
    path_to_tsv: str
    slice_idx : str, optional
        Index level value to slice on, if required to remove potentially function-breaking data. Requires `slice_lvl` to also be defined.
    slice_lvl : str, optional
        Index level name to slice on, if required to remove potentially function-breaking data. Requires `slice_idx` to also be defined.
    """
    df = pd.read_csv(path_to_tsv, delimiter='\t', index_col=0)
    index_names = df.index.name.split("\\")[0].split(",")
    column_name = df.index.name.split("\\")[1]
    df.index = df.index.str.split(',', expand=True).rename(index_names)
    df.columns = df.columns.rename(column_name)
    if df.columns.name == "time":
        df.columns = df.columns.astype(int)
    else:
        df.columns = df.columns.str.strip()

    if slice_idx is not None:
        df = df.xs(slice_idx, level=slice_lvl)

    return df.apply(to_numeric)


def remove_digits():
    """
    Functionality to be passed to str.translate to remove numbers from
    string endings
    """
    return str.maketrans("", "", digits)


def fill_spatio_temporal_data_gaps(
    data_with_gaps: xr.DataArray, gap_filling_methods: dict, data_is_relative: bool = True
) -> xr.DataArray:


    def _temporal_gap_filling(data):
        if gap_filling_methods["year-gap-filling-method"] == "average":
            missing_year_filler = data.mean("year")
            if data_is_relative:
                missing_year_filler = missing_year_filler / missing_year_filler.sum("end_use")
            data = data.fillna(missing_year_filler)
        else:
            interpolation_method = gap_filling_methods["year-gap-filling-method"]
            data = data.interpolate_na(
                "year", method=interpolation_method, fill_value="extrapolate"
            )
        return data

    def _spatial_gap_filling(data):
        # ASSUME: missing countries have the same percentage contribution as their neighbours
        for country_group_mappings in gap_filling_methods["country-gap-filling"].values():
            country_averages = (
                data
                .sel(country_code=country_group_mappings["countries_with_data"])
                .mean("country_code")
            )
            if data_is_relative:
                country_averages = country_averages / country_averages.sum("end_use")
            for country in country_group_mappings["countries_without_data"]:
                if country not in data.country_code:
                    data = merge_da(
                        [data, country_averages.expand_dims(country_code=[country])]
                    )
        return data

    data_without_temporal_gaps = _temporal_gap_filling(data_with_gaps)
    data_without_gaps = _spatial_gap_filling(data_without_temporal_gaps)

    # Sometimes adding in new countries leads to new data gaps be be filled in
    if data_without_gaps.isnull().any():
        data_without_temporal_gaps = _temporal_gap_filling(data_without_gaps)
        data_without_gaps = _temporal_gap_filling(data_without_temporal_gaps)

    return data_without_gaps
