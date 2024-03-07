"""Utility functions."""
import pycountry

from string import digits


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(eu_country_code) == 2, "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)

    return convert_country_code(eu_country_code, output="alpha3")


def convert_country_code(input_country, output="alpha3"):
    """
    Converts input country code or name into either a 2- or 3-letter code.

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


# conversion utils
def ktoe_to_twh(array):
    """Convert KTOE to TWH"""
    return array * 1.163e-2


def pj_to_twh(array):
    """Convert PJ to TWh"""
    return array / 3.6


def tj_to_twh(array):
    """Convert TJ to TWh"""
    return pj_to_twh(array) / 1000


def gwh_to_tj(array):
    """Convert GWh to TJ"""
    return array * 3.6


def remove_digits():
    """
    Functionality to be passed to str.translate to remove numbers from
    string endings
    """
    return str.maketrans("", "", digits)
