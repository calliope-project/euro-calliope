"""Utility functions."""
import pycountry
import pandas as pd


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert (
        len(eu_country_code) == 2
    ), "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)
    if eu_country_code.lower() == "el":
        iso2 = "gr"
    elif eu_country_code.lower() == "uk":
        iso2 = "gb"
    elif (
        eu_country_code.lower() == "bh"
    ):  # this is a weird country code used in the biofuels dataset
        iso2 = "ba"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3


def get_alpha2(country, eurostat=True):
    """
    Returns the alpha-2 ISO country code for the given country.
    """
    if country in ["United Kingdom", "GB", "GBR"] and eurostat is True:
        return "UK"
    elif country in ["Greece", "GR", "GRC"] and eurostat is True:
        return "EL"
    else:
        return pycountry.countries.lookup(country).alpha_2


def to_numeric(series):
    """
    Clean up a pandas.Series which was parsed as strings, but is really numeric:

    1. replace "-" for "NaN" into numbers and NaNs
    2. removes random superscript attached to numbers
       (e.g. pointing to footnotes in an excel), "1000c" -> 1000

    Returns a numeric pandas.Series.

    """
    series = series.astype(str).str.extract("(\-*\d+\.*\d*)")[0]
    return pd.to_numeric(series, errors="coerce")


def gwh_to_tj(array):
    """Convert GWh to TJ"""
    return array * 3.6

def pj_to_twh(array):
    """Convert PJ to TWh"""
    return array / 3.6

def tj_to_twh(array):
    """Convert TJ to TWh"""
    return pj_to_twh(array) * 1e-3

def read_tdf(filename):
    """ 
    Read a tidy dataframe from CSV. 
    This assumes that all except the final column is an index level,  so returns a MultiIndexed Pandas Series. 
    """
    df = pd.read_csv(filename, header=0)
    tdf = df.set_index([i for i in df.columns[:-1]]).squeeze()
    return tdf
