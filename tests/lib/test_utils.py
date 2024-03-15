import pytest
from eurocalliopelib.utils import convert_country_code, eu_country_code_to_iso3


@pytest.mark.parametrize(
    "eu_country_code,iso3166", [("uk", "GBR"), ("el", "GRC"), ("de", "DEU")]
)
def test_eu_country_code_conversion(eu_country_code, iso3166):
    assert eu_country_code_to_iso3(eu_country_code) == iso3166


@pytest.mark.parametrize(
    "input,alpha3",
    [
        ("uk", "GBR"),
        ("el", "GRC"),
        ("de", "DEU"),
        ("Germany", "DEU"),
        ("Greece", "GRC"),
    ],
)
def test_conversion_to_alpha3(input, alpha3):
    assert convert_country_code(input, "alpha3") == alpha3


@pytest.mark.parametrize(
    "input,alpha2",
    [("uk", "GB"), ("el", "GR"), ("de", "DE"), ("Germany", "DE"), ("Greece", "GR")],
)
def test_conversion_to_alpha2(input, alpha2):
    assert convert_country_code(input, "alpha2") == alpha2


@pytest.mark.parametrize(
    "input,alpha2",
    [("uk", "UK"), ("el", "EL"), ("de", "DE"), ("Germany", "DE"), ("Greece", "EL")],
)
def test_conversion_to_eu_alpha2(input, alpha2):
    assert convert_country_code(input, "alpha2_eu") == alpha2
