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
    elif eu_country_code.lower() == "bh": # this is a weird country code used in the biofuels dataset
        iso2 = "ba"
    else:
        iso2 = eu_country_code
    return pycountry.countries.lookup(iso2).alpha_3
