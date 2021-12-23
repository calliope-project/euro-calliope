import pytest

from eurocalliopelib import filters


@pytest.mark.parametrize('value,unit_string,expected_string', [
    (4.0, "MWh", "(4 MWh)"),
    (4.2, "MWh", "(4.2 MWh)"),
    (4.44445, "MWh", "(4.44445 MWh)"),
    (3, "MW", "(3 MW)"),
    (10000, "W", "(10,000 W)"),
    (100000, "EUR2015/MW/year", "(100,000 EUR2015/MW/year)")
])
def test_unit_filter(value, unit_string, expected_string):
    assert filters.unit(value, unit_string) == expected_string


def test_unit_filter_without_parenthesis():
    assert filters.unit(4.0, "MWh", parenthesis=False) == "4 MWh"
