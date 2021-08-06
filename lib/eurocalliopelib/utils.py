"""Utility functions."""
import pycountry
import pandas as pd
from datetime import datetime

from string import digits

UNIT_CONVERSION_MAPPING = {
    ("gwh", "tj"): 3.6,
    ("ktoe", "twh"): 1.163e-2,
    ("tj", "twh"): 1e3 / 3.6,
    ("pj", "twh"): 1 / 3.6,
    ("tj", "ktoe"): 23.88e-3
}


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


def to_numeric(series):
    """
    Clean up a pandas.Series which was parsed as strings, but is really numeric:

    1. replace "-" for "NaN" into numbers and NaNs
    2. removes random superscript attached to numbers
       (e.g. pointing to footnotes in an excel), "1000c" -> 1000

    Returns a numeric pandas.Series.

    """
    series_name = series.name
    series = series.astype(str).str.extract("(\\-*\\d+\\.*\\d*)")[0]
    return pd.to_numeric(series.rename(series_name), errors="coerce")


def convert_unit(df, output_unit, input_unit=None, unit_in_output_idx=True):
    df = df.copy()
    if isinstance(df.index, pd.MultiIndex) and "unit" in df.index.names:
        units = df.index.get_level_values("unit").unique()
        if input_unit is None:
            assert len(units) == 1, f"Cannot infer unit for data with multiple available units {units}"
            input_unit = units[0]
        mask = df.index.get_level_values("unit") == input_unit
    else:
        units = None
        if input_unit is None:
            raise ValueError("Cannot infer unit for data")
        else:
            mask = df.index

    if input_unit.lower() != output_unit.lower():
        df.loc[mask] *= UNIT_CONVERSION_MAPPING[(input_unit.lower(), output_unit.lower())]

    if units is None and unit_in_output_idx is True:
        df = add_idx_level(df, unit=output_unit.lower())
    elif units is not None:
        if unit_in_output_idx is True:
            idx_renamer = lambda x: output_unit.lower() if x == input_unit else x
            df = df.rename(idx_renamer, level="unit")
        else:
            assert len(units) == 1, f"Cannot drop the index level `unit` with multiple units {units}"
            df = df.xs(input_unit, level="unit")
    return df


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
        please_squeeze = True
    else:
        please_squeeze = False
        orig_cols = new_data.columns
    new_data = (
        new_data
        .assign(**level_info)
        .set_index([i for i in level_info.keys()], append=True)
    )
    if please_squeeze:
        new_data = new_data.squeeze(axis=1)
    else:
        new_data.columns = orig_cols
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
