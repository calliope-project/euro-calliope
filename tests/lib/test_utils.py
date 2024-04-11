import numpy as np
import pandas as pd
import pytest
import xarray as xr
from eurocalliopelib.utils import (
    convert_country_code,
    convert_valid_countries,
    rename_and_groupby,
)


class TestRenameAndGroupby:
    @pytest.fixture
    def da(self, request):
        data = [1, 2, 3, 4, 5]
        if request.param == "single_index":
            index = pd.Index(["A", "B", "C", "D", "E"], name="foo")
        elif request.param == "multi_index":
            data += data
            index = pd.MultiIndex.from_product(
                [["X", "Y"], ["A", "B", "C", "D", "E"]], names=["bar", "foo"]
            )
        return pd.Series(data, index=index).to_xarray()

    @pytest.fixture
    def da_with_nan(self):
        data = [1, np.nan, np.nan, 4, 5]
        data += reversed(data)
        index = pd.MultiIndex.from_product(
            [["X", "Y"], ["A", "B", "C", "D", "E"]], names=["bar", "foo"]
        )
        return pd.Series(data, index=index).to_xarray()

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"A": "A1", "B": "B1", "C": "C1", "D": "D1", "E": "E1"},
                {"A1": 1, "B1": 2, "C1": 3, "D1": 4, "E1": 5},
            ),
            ({"A": "A1", "B": "A1", "C": "A1", "D": "A1", "E": "A1"}, {"A1": 15}),
            (
                {"A": "A1", "B": "A1", "C": "B1", "D": "B1", "E": "E1"},
                {"A1": 3, "B1": 7, "E1": 5},
            ),
            ({"A": "A1", "B": "B1"}, {"A1": 1, "B1": 2}),
            ({"A": "A1", "B": "B1", "C": "B1"}, {"A1": 1, "B1": 5}),
        ],
    )
    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_groupby(self, rename_dict, expected, da):
        new_da = rename_and_groupby(
            da, rename_dict, keep_non_renamed=False, dropna=False, dim_name="foo"
        )

        assert (
            new_da.coords["foo"].to_index().symmetric_difference(expected.keys()).empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(foo=dim_item), da_data)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"A": "A1", "B": "B1", "C": "C1", "D": "D1", "E": "E1"},
                {"A1": 1, "B1": 2, "C1": 3, "D1": 4, "E1": 5},
            ),
            ({"A": "A1", "B": "A1", "C": "A1", "D": "A1", "E": "A1"}, {"A1": 15}),
            (
                {"A": "A1", "B": "A1", "C": "B1", "D": "B1", "E": "E1"},
                {"A1": 3, "B1": 7, "E1": 5},
            ),
            ({"A": "A1", "B": "B1"}, {"A1": 1, "B1": 2, "C": 3, "D": 4, "E": 5}),
            ({"A": "A1", "B": "B1", "C": "B1"}, {"A1": 1, "B1": 5, "D": 4, "E": 5}),
        ],
    )
    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_groupby_keep_renamed(self, rename_dict, expected, da):
        new_da = rename_and_groupby(
            da, rename_dict, keep_non_renamed=True, dim_name="foo"
        )

        assert (
            new_da.coords["foo"].to_index().symmetric_difference(expected.keys()).empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(foo=dim_item), da_data)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"A": "A1", "B": "B1", "C": "C1", "D": "D1", "E": "E1"},
                {"A1": [1, 5], "B1": [np.nan, 4], "D1": [4, np.nan], "E1": [5, 1]},
            ),
            ({"A": "A1", "B": "A1", "C": "A1", "D": "A1", "E": "A1"}, {"A1": [10, 10]}),
            (
                {"A": "A1", "B": "A1", "C": "B1", "D": "B1", "E": "E1"},
                {"A1": [1, 9], "B1": [4, np.nan], "E1": [5, 1]},
            ),
            ({"A": "A1", "B": "B1"}, {"A1": [1, 5], "B1": [np.nan, 4]}),
            ({"A": "A1", "B": "B1", "C": "B1"}, {"A1": [1, 5], "B1": [np.nan, 4]}),
        ],
    )
    def test_rename_groupby_dropna(self, da_with_nan, rename_dict, expected):
        new_da = rename_and_groupby(
            da_with_nan, rename_dict, dim_name="foo", dropna=True
        )

        assert (
            new_da.coords["foo"].to_index().symmetric_difference(expected.keys()).empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(foo=dim_item), da_data)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"A": "A1", "B": "B1", "C": "C1", "D": "D1", "E": "E1"},
                {
                    "A1": [1, 5],
                    "B1": [np.nan, 4],
                    "C1": [np.nan, np.nan],
                    "D1": [4, np.nan],
                    "E1": [5, 1],
                },
            ),
            ({"A": "A1", "B": "A1", "C": "A1", "D": "A1", "E": "A1"}, {"A1": [10, 10]}),
            (
                {"A": "A1", "B": "A1", "C": "B1", "D": "B1", "E": "E1"},
                {"A1": [1, 9], "B1": [4, np.nan], "E1": [5, 1]},
            ),
            ({"A": "A1", "B": "B1"}, {"A1": [1, 5], "B1": [np.nan, 4]}),
            ({"A": "A1", "B": "B1", "C": "B1"}, {"A1": [1, 5], "B1": [np.nan, 4]}),
        ],
    )
    def test_rename_groupby_no_dropna(self, da_with_nan, rename_dict, expected):
        new_da = rename_and_groupby(
            da_with_nan, rename_dict, dim_name="foo", dropna=False
        )

        assert (
            new_da.coords["foo"].to_index().symmetric_difference(expected.keys()).empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(foo=dim_item), da_data)

    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_dim(self, da):
        rename_dict = {"A": "A1", "B": "B1", "C": "C1", "D": "D1", "E": "E1"}
        expected = {"A1": 1, "B1": 2, "C1": 3, "D1": 4, "E1": 5}
        new_da = rename_and_groupby(
            da, rename_dict, dim_name="foo", new_dim_name="baz", dropna=False
        )

        assert "baz" in new_da.coords
        assert "foo" not in new_da.coords
        assert (
            new_da.coords["baz"].to_index().symmetric_difference(expected.keys()).empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(baz=dim_item), da_data)


@pytest.mark.parametrize(
    "input_country,iso3166,output",
    [
        ("uk", "GBR", "alpha3"),
        ("el", "GRC", "alpha3"),
        ("de", "DEU", "alpha3"),
        ("Germany", "DE", "alpha2"),
        ("United Kingdom", "GB", "alpha2"),
        ("United Kingdom", "UK", "alpha2_eurostat"),
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
        ("United Kingdom", "GBR"),
    ],
)
def test_country_code_conversion_no_output_specified(input_country, iso3166):
    assert convert_country_code(input_country) == iso3166


@pytest.mark.parametrize(
    ["countries", "expected_output"],
    [
        (["uk", "fr"], {"uk": "GBR", "fr": "FRA"}),
        (["el", "ba"], {"el": "GRC", "ba": "BIH"}),
        (["foo", "bar"], {}),
        (["Germany", "foo"], {"Germany": "DEU"}),
        (
            xr.DataArray(pd.Index(["Germany", "foo"], name="bar")).bar.values,
            {"Germany": "DEU"},
        ),
    ],
)
def test_valid_list_of_country_codes(countries, expected_output):
    mapped_countries = convert_valid_countries(countries)
    assert mapped_countries == expected_output


@pytest.mark.parametrize(
    "input,alpha2",
    [("uk", "UK"), ("el", "EL"), ("de", "DE"), ("Germany", "DE"), ("Greece", "EL")],
)
def test_conversion_to_eu_alpha2(input, alpha2):
    assert convert_country_code(input, "alpha2_eurostat") == alpha2
