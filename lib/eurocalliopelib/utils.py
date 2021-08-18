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


def rename_and_groupby(da, rename_dict, dim, new_group_dim_name=None, dropna=False):
    if new_group_dim_name is None:
        new_group_dim_name = dim
    rename_da = xr.DataArray(pd.Series(rename_dict).rename_axis(index=dim).rename(new_group_dim_name))
    da = (
        da
        .reindex({dim: rename_da[dim]})
        .groupby(rename_da)
        .sum(new_group_dim_name, skipna=True, min_count=1, keep_attrs=True)
    )
    if dropna:
        da = da.dropna(new_group_dim_name, how="all")
    return da


def merge_da(da_list):
    datasets = [da.to_dataset(name="var") for da in da_list]
    return xr.merge(datasets, combine_attrs="no_conflicts")["var"]


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
    """Convert GWh to TWh"""
    return array / 3.6


def tj_to_twh(array):
    """Convert GWh to TWh"""
    return pj_to_twh(array) / 1000


def ktoe_to_twh(array):
    """Convert KTOE to TWH"""
    return array * 1.163e-2


def add_idx_level(data, **level_info):
    """
    Add a new index level or multiple levels to a Pandas dataset (Series or DataFrame).

    Parameters
    ----------
    data : pd.Series or pd.DataFrame
        data to which Index level is added
    level_info :
        kwargs where keys are the level names, and values are the level values.
    """
    new_data = data.copy()
    if isinstance(new_data, pd.Series):
        if new_data.name in level_info.keys():
            raise KeyError(
                "Cannot add index level of the same name as the pandas Series"
            )
        new_data = new_data.to_frame()
    new_data = (
        new_data
        .assign(**level_info)
        .set_index([i for i in level_info.keys()], append=True)
        .squeeze(axis=1)
    )

    return new_data


def read_eurostat_tsv(path_to_tsv, index_names, slice_idx=None, slice_lvl=None):
    """

    Read a typical tab-delimited file from EUROSTAT. These have a specific structure
    where the data is tab-delimited but the multi-index data is comma delimited.
    This function also prepares the data in the expectation that it is all numeric
    and that the columns are given as years (a standard EUROSTAT format)

    Parameters
    ---------
    path_to_tsv: str
    index_names : array
        names to which the index levels correspond. Must be the same length as
        the number of expected index levels.
    slice_idx : str, optional
        Index level value to slice on, if required to remove potentially function-breaking data. Requires `slice_lvl` to also be defined.
    slice_lvl : str, optional
        Index level name to slice on, if required to remove potentially function-breaking data. Requires `slice_idx` to also be defined.
    """
    df = pd.read_csv(path_to_tsv, delimiter='\t', index_col=0)
    df.index = df.index.str.split(',', expand=True).rename(index_names)
    if slice_idx is not None:
        df = df.xs(slice_idx, level=slice_lvl)
    df.columns = df.columns.astype(int).rename("year")
    return df.apply(to_numeric)


def remove_digits():
    """
    Functionality to be passed to str.translate to remove numbers from
    string endings
    """
    return str.maketrans("", "", digits)
