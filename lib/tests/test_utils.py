import pytest

from eurocalliopelib.utils import eu_country_code_to_iso3


@pytest.mark.parametrize('eu_country_code,iso3166', [
    ("uk", "GBR"),
    ("el", "GRC"),
    ("de", "DEU")
])
def test_country_code_conversion(eu_country_code, iso3166):
    assert eu_country_code_to_iso3(eu_country_code) == iso3166
