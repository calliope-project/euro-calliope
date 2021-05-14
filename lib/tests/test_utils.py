import pytest

from eurocalliopelib.utils import convert_country_code


@pytest.mark.parametrize(
    "input_country,iso3166,output",
    [
        ("uk", "GBR", "alpha3"),
        ("el", "GRC", "alpha3"),
        ("de", "DEU", "alpha3"),
        ("Germany", "DE", "alpha2"),
        ("Great Britain", "GB", "alpha2"),
        ("Great Britain", "UK", "alpha2_eurostat"),
        ("Greece", "GR", "alpha2"),
        ("Greece", "EL", "alpha2_eurostat"),
        ("bh", "BA", "alpha2"),
    ],
)
def test_country_code_conversion(input_country, iso3166, output):
    assert convert_country_code(input_country, output) == iso3166


@pytest.mark.parametrize(
    "input_country,iso3166",
    [
        ("uk", "GBR"),
        ("el", "GRC"),
        ("de", "DEU"),
        ("Germany", "DEU"),
        ("Great Britain", "GBR"),
    ],
)
def test_country_code_conversion_no_output_specified(input_country, iso3166):
    assert convert_country_code(input_country) == iso3166
