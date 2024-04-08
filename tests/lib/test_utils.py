import numpy as np
import pandas as pd
import pytest
from eurocalliopelib.utils import (
    convert_country_code,
    merge_da,
    remove_digits,
    rename_and_groupby,
    to_numeric,
)
from pandas.api.types import is_numeric_dtype


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


class TestMergeDA:
    @pytest.fixture
    def da1(self):
        return pd.Series(
            [1, 2, 3, 4],
            index=pd.MultiIndex.from_product(
                [["X", "Y"], ["A", "B"]], names=["bar", "foo"]
            ),
        ).to_xarray()

    @pytest.fixture
    def da2(self):
        return pd.Series(
            [1], index=pd.MultiIndex.from_product([["Z"], ["C"]], names=["bar", "foo"])
        ).to_xarray()

    @pytest.fixture
    def da3(self):
        return pd.Series(
            [1, 2],
            index=pd.MultiIndex.from_product([["X", "Y"], ["D"]], names=["bar", "foo"]),
        ).to_xarray()

    @pytest.fixture
    def da4(self):
        return pd.Series(
            [10, 20, 30, 40],
            index=pd.MultiIndex.from_product(
                [["X", "Y"], ["C"], ["P", "Q"]], names=["bar", "foo", "foobar"]
            ),
        ).to_xarray()

    @staticmethod
    def assertions(new_da, old_das, bar_dims, foo_dims, data_sum):
        assert new_da["bar"].to_index().symmetric_difference(bar_dims).empty
        assert new_da["foo"].to_index().symmetric_difference(foo_dims).empty
        for old_da in old_das:
            np.testing.assert_array_equal(new_da.sel(old_da.coords), old_da)
        assert new_da.sum() == data_sum

    def test_merge_da_extend_two_dims(self, da1, da2):
        new_da = merge_da([da1, da2])
        self.assertions(new_da, [da1, da2], ["X", "Y", "Z"], ["A", "B", "C"], 11)

    def test_merge_da_extend_one_dim(self, da1, da3):
        new_da = merge_da([da1, da3])
        self.assertions(new_da, [da1, da3], ["X", "Y"], ["A", "B", "D"], 13)

    def test_merge_da_extend_one_dim_and_hanging_dim(self, da1, da4):
        # for arrays without a particular dimension, data on other dimensions are repeated on the new dimension.
        # In this case, we get da1 data going from:
        #   bar  foo
        #   X    A      1
        #        B      2
        #   Y    A      3
        #        B      4
        # to:
        #   bar  foo  foobar
        #   X    A    P         1.0
        #             Q         1.0
        #        B    P         2.0
        #             Q         2.0
        #   Y    A    P         3.0
        #             Q         3.0
        #        B    P         4.0
        #             Q         4.0

        new_da = merge_da([da1, da4])
        self.assertions(new_da, [da4], ["X", "Y"], ["A", "B", "C"], 120)

        assert new_da["foobar"].to_index().symmetric_difference(["P", "Q"]).empty
        np.testing.assert_array_equal(new_da.sel(da4.coords), da4)

    def test_merge_da_extend_two_dims_with_three_das(self, da1, da2, da3):
        new_da = merge_da([da1, da2, da3])
        self.assertions(
            new_da, [da1, da2, da3], ["X", "Y", "Z"], ["A", "B", "C", "D"], 14
        )

    @pytest.mark.parametrize("new_name", [None, "foo"])
    def test_rename_merge_da(self, da1, da2, new_name):
        new_da = merge_da([da1, da2], merged_da_name=new_name)
        if new_name is None:
            assert new_da.name is None
        else:
            assert new_da.name == new_name


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


@pytest.mark.parametrize("string_with_digits", ["trial01"])
def test_digit_removal_from_string(string_with_digits):
    assert string_with_digits.translate(remove_digits()) == "trial"


class TestToNumeric:
    def test_dash_to_NaN(self):
        new_series = to_numeric(pd.Series(["-", 1, "-"]))
        assert new_series.isna().sum() == 2
        assert is_numeric_dtype(new_series)

    @pytest.mark.parametrize(
        ("raw", "expected"),
        [
            ([1, "1a", "2-"], [1, 1, 2]),
            ([10, "100a", "1000-"], [10, 100, 1000]),
            ([1.0, "1.0a", "-1-"], [1.0, 1.0, -1]),
        ],
    )
    def test_remove_non_numeric_characters(self, raw, expected):
        new_series = to_numeric(pd.Series(raw))
        assert np.allclose(new_series, expected)
        assert is_numeric_dtype(new_series)

    def test_no_change_to_numeric_characters(self):
        original_series = pd.Series([0, 1, 100, 2.0, np.nan])
        new_series = to_numeric(original_series.copy())
        assert new_series.equals(original_series)
