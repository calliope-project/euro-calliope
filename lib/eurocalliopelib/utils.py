"""Utility functions."""
import pycountry
import pandas as pd


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
