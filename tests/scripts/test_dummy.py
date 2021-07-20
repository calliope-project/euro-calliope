# This is a proof of concept. It should be removed once other tests are added.
from scripts.biofuels import utils # This proofs script imports are possible.


def test_country_code_conversion():
    assert utils.eu_country_code_to_iso3("uk") == "GBR"
