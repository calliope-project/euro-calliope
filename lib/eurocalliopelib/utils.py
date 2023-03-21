"""Utility functions."""
import pycountry


def eu_country_code_to_iso3(eu_country_code):
    """Converts EU country code to ISO 3166 alpha 3.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(eu_country_code) == 2, "EU country codes are of length 2, yours is '{}'.".format(eu_country_code)
    if eu_country_code.lower() == "el":
        iso2 = "gr"
    elif eu_country_code.lower() == "uk":
        iso2 = "gb"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3


def iso3_to_eu_country_code(iso3_code):
    """Converts ISO 3166 alpha 3 to EU country code.
    The European Union uses its own country codes, which often but not always match ISO 3166.
    """
    assert len(iso3_code) == 3, "ISO 3166 alpha 3 country codes are of length 3, yours is '{}'.".format(iso3_code)
    if iso3_code.lower() == "grc":
        eu_country_code = "el"
    elif iso3_code.lower() == "gbr":
        eu_country_code = "uk"
    else:
        eu_country_code = pycountry.countries.lookup(iso3_code).alpha_2
    return eu_country_code
