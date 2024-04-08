"""Utility functions."""

from string import digits
from typing import Optional

# We import netCDF4 before xarray to mitigate a numpy warning that translates to an error on linux:
# https://github.com/pydata/xarray/issues/7259
import netCDF4  # noqa: F401
import pandas as pd
import pycountry
import xarray as xr


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert (
        len(eu_country_code) == 2
    ), f"EU country codes are of length 2, yours is '{eu_country_code}'."

    return convert_country_code(eu_country_code, output="alpha3")


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


def rename_and_groupby(
    da: xr.DataArray,
    rename_dict: dict,
    dim_name: str,
    new_dim_name: Optional[str] = None,
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
        new_dim_name (Optional[str], optional): Defaults to None.
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
        rename_series = rename_series.reindex(existing_dim_items).fillna(
            existing_dim_items
        )

    if new_dim_name is None:
        new_dim_name = f"_{dim_name}"  # placeholder that we'll revert
        revert_dim_name = True
    else:
        revert_dim_name = False

    rename_da = xr.DataArray(rename_series.rename(new_dim_name))
    da = (
        da.reindex({dim_name: rename_da[dim_name]})
        .groupby(rename_da)
        .sum(dim_name, skipna=True, min_count=1, keep_attrs=True)
    )
    if revert_dim_name:
        da = da.rename({new_dim_name: dim_name})
        new_dim_name = dim_name
    if dropna:
        da = da.dropna(new_dim_name, how="all")
    return da


def merge_da(da_list: list, merged_da_name: Optional[str] = None) -> xr.DataArray:
    """
    Merge dataArrays with the same dimensions but different dimension items
    into a single xarray datarray

    Args:
        da_list (list): List of xarray DataArrays.
        merged_da_name (Optional[str], optional): Name of merged datarray. Defaults to None.

    Returns:
        xr.DataArray:
            Merged Datarray, in which all dimensions contain all items defined in the
            arrays in `da_list`

    """
    datasets = [da.rename("var") for da in da_list]
    return xr.merge(datasets, combine_attrs="no_conflicts")["var"].rename(
        merged_da_name
    )


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


def pj_to_twh(array):
    """Convert PJ to TWh"""
    return array / 3.6


def tj_to_twh(array):
    """Convert TJ to TWh"""
    return pj_to_twh(array) / 1000


def ktoe_to_twh(array):
    """Convert KTOE to TWH"""
    return array * 1.163e-2


def gwh_to_tj(array):
    """Convert GWh to TJ"""
    return array * 3.6


def remove_digits():
    """
    Functionality to be passed to str.translate to remove numbers from
    string endings
    """
    return str.maketrans("", "", digits)
