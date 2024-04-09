"""Utility functions."""

import pycountry


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
    Converts input country code or name into either a 2- or 3-letter code.

    ISO alpha2: alpha2
    ISO alpha2 with Eurostat codes: alpha2_eu
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

    if output == "alpha2_eu":
        result = pycountry.countries.lookup(input_country).alpha_2
        if result == "GB":
            return "UK"
        elif result == "GR":
            return "EL"
        else:
            return result

    if output == "alpha3":
        return pycountry.countries.lookup(input_country).alpha_3


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


def tj_to_ktoe(array):
    """Convert TJ to Ktoe"""
    return array * 23.88e-3
