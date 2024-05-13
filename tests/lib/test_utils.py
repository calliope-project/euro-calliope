import logging

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
        if request.param == "single_index":
            data = [1, 2, 3, 4, 5]
            index = pd.Index(["FRA", "DEU", "GBR", "ITA", "ESP"], name="country_code")
        elif request.param == "multi_index":
            # Same values on both "other_dim" levels, which allows us to use
            # `np.testing.assert_array_equal` for both `multi_index` and `single_index`
            # because both `np.testing.assert_array_equal([1], 1)` and `np.testing.assert_array_equal([1, 1], 1)`
            # are valid matches.
            data = [1, 2, 3, 4, 5, 1, 2, 3, 4, 5]
            index = pd.MultiIndex.from_product(
                [["X", "Y"], ["FRA", "DEU", "GBR", "ITA", "ESP"]],
                names=["other_dim", "country_code"],
            )
        return pd.Series(data, index=index).to_xarray()

    @pytest.fixture
    def da_with_nan(self):
        # Only country with NaN in both "other_dim" values is GBR.
        data = [1, np.nan, np.nan, 4, 5, 5, 4, np.nan, np.nan, 1]
        index = pd.MultiIndex.from_product(
            [["X", "Y"], ["FRA", "DEU", "GBR", "ITA", "ESP"]],
            names=["other_dim", "country_code"],
        )
        return pd.Series(data, index=index).to_xarray()

    @staticmethod
    def assert_expected_dims_and_dim_vals(new_da: xr.DataArray, expected: dict):
        # check that the dimension items have been appropriately renamed
        assert (
            new_da.coords["country_code"]
            .to_index()
            .symmetric_difference(expected.keys())
            .empty
        )
        # check that the values associated with each dimension item have been appropriated grouped
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(country_code=dim_item), da_data)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"FRA": "FR", "DEU": "DE", "GBR": "GB", "ITA": "IT", "ESP": "ES"},
                {"FR": 1, "DE": 2, "GB": 3, "IT": 4, "ES": 5},
            ),
            (
                {"FRA": "FR", "DEU": "FR", "GBR": "FR", "ITA": "FR", "ESP": "FR"},
                {"FR": 15},
            ),
            (
                {"FRA": "FR", "DEU": "FR", "GBR": "DE", "ITA": "DE", "ESP": "ES"},
                {"FR": 3, "DE": 7, "ES": 5},
            ),
        ],
    )
    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    @pytest.mark.parametrize("drop_other_dim_items", [True, False])
    def test_rename_groupby(self, rename_dict, expected, da, drop_other_dim_items):
        "`drop_other_dim_items` doesn't do anything in these examples as all dimension items are explicitly referenced"
        new_da = rename_and_groupby(
            da,
            rename_dict,
            drop_other_dim_items=drop_other_dim_items,
            dropna=False,
            dim_name="country_code",
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            ({"FRA": "FR", "DEU": "DE"}, {"FR": 1, "DE": 2}),
            ({"FRA": "FR", "DEU": "DE", "GBR": "DE"}, {"FR": 1, "DE": 5}),
        ],
    )
    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_groupby_drop_other_dim_items(self, rename_dict, expected, da):
        new_da = rename_and_groupby(
            da,
            rename_dict,
            drop_other_dim_items=True,
            dropna=False,
            dim_name="country_code",
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"FRA": "FR", "DEU": "DE"},
                {"FR": 1, "DE": 2, "GBR": 3, "ITA": 4, "ESP": 5},
            ),
            (
                {"FRA": "FR", "DEU": "DE", "GBR": "DE"},
                {"FR": 1, "DE": 5, "ITA": 4, "ESP": 5},
            ),
        ],
    )
    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_groupby_keep_other_dim_items(self, rename_dict, expected, da):
        new_da = rename_and_groupby(
            da, rename_dict, drop_other_dim_items=False, dim_name="country_code"
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"FRA": "FR", "DEU": "FR", "GBR": "FR", "ITA": "FR", "ESP": "FR"},
                # Multiple values expected now as we have different values on the "other_dim" dimension.
                {"FR": [10, 10]},
            ),
            (
                {"FRA": "FR", "DEU": "FR", "GBR": "DE", "ITA": "DE", "ESP": "ES"},
                {"FR": [1, 9], "DE": [4, np.nan], "ES": [5, 1]},
            ),
            ({"FRA": "FR", "DEU": "DE"}, {"FR": [1, 5], "DE": [np.nan, 4]}),
            (
                {"FRA": "FR", "DEU": "DE", "GBR": "DE"},
                {"FR": [1, 5], "DE": [np.nan, 4]},
            ),
        ],
    )
    @pytest.mark.parametrize("dropna", [True, False])
    def test_rename_groupby_with_na(self, da_with_nan, rename_dict, expected, dropna):
        "dropna doesn't do anything in these examples as `how=all` in the underlying call"
        new_da = rename_and_groupby(
            da_with_nan, rename_dict, dim_name="country_code", dropna=dropna
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"FRA": "FR", "DEU": "DE", "GBR": "GB", "ITA": "IT", "ESP": "ES"},
                # GB dimension dropped
                {"FR": [1, 5], "DE": [np.nan, 4], "IT": [4, np.nan], "ES": [5, 1]},
            ),
        ],
    )
    def test_rename_groupby_dropna_removes_dim(
        self, da_with_nan, rename_dict, expected
    ):
        new_da = rename_and_groupby(
            da_with_nan, rename_dict, dim_name="country_code", dropna=True
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize(
        ("rename_dict", "expected"),
        [
            (
                {"FRA": "FR", "DEU": "DE", "GBR": "GB", "ITA": "IT", "ESP": "ES"},
                {
                    "FR": [1, 5],
                    "DE": [np.nan, 4],
                    "GB": [np.nan, np.nan],  # this dimension stays
                    "IT": [4, np.nan],
                    "ES": [5, 1],
                },
            ),
        ],
    )
    def test_rename_groupby_no_dropna_keeps_dim(
        self, da_with_nan, rename_dict, expected
    ):
        new_da = rename_and_groupby(
            da_with_nan, rename_dict, dim_name="country_code", dropna=False
        )
        self.assert_expected_dims_and_dim_vals(new_da, expected)

    @pytest.mark.parametrize("da", ["single_index", "multi_index"], indirect=True)
    def test_rename_dim(self, da):
        rename_dict = {"FRA": "FR", "DEU": "DE", "GBR": "GB", "ITA": "IT", "ESP": "ES"}
        expected = {"FR": 1, "DE": 2, "GB": 3, "IT": 4, "ES": 5}
        new_da = rename_and_groupby(
            da,
            rename_dict,
            dim_name="country_code",
            new_dim_name="country",
            dropna=False,
        )

        assert "country" in new_da.coords
        assert "country_code" not in new_da.coords
        assert (
            new_da.coords["country"]
            .to_index()
            .symmetric_difference(expected.keys())
            .empty
        )
        for dim_item, da_data in expected.items():
            np.testing.assert_array_equal(new_da.sel(country=dim_item), da_data)


@pytest.mark.parametrize(
    "input_country,iso3166,output",
    [
        ("uk", "GBR", "alpha3"),
        ("el", "GRC", "alpha3"),
        ("de", "DEU", "alpha3"),
        ("Germany", "DE", "alpha2"),
        ("United Kingdom", "GB", "alpha2"),
        ("United Kingdom", "UK", "alpha2_eu"),
        ("Greece", "GR", "alpha2"),
        ("Greece", "EL", "alpha2_eu"),
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
    ["countries", "expected_output", "expected_log"],
    [
        (["uk", "fr"], {"uk": "GBR", "fr": "FRA"}, False),
        (["el", "ba"], {"el": "GRC", "ba": "BIH"}, False),
        (["foo", "bar"], {}, True),
        (["Germany", "foo"], {"Germany": "DEU"}, True),
        (
            xr.DataArray(pd.Index(["Germany", "foo"], name="bar")).bar.values,
            {"Germany": "DEU"},
            True,
        ),
    ],
)
def test_valid_list_of_country_codes(caplog, countries, expected_output, expected_log):
    caplog.set_level(logging.INFO)
    mapped_countries = convert_valid_countries(countries, errors="ignore")
    assert mapped_countries == expected_output
    assert ("Skipping country/region" in caplog.text) is expected_log


@pytest.mark.parametrize(
    ["countries", "expected_error"],
    [
        (["uk", "fr"], False),
        (["el", "ba"], False),
        (["foo", "bar"], True),
        (["Germany", "foo"], True),
    ],
)
def test_valid_list_of_country_codes_catch_errors(countries, expected_error):
    if expected_error:
        with pytest.raises(LookupError):
            convert_valid_countries(countries, errors="raise")
    else:
        convert_valid_countries(countries, errors="raise")


@pytest.mark.parametrize(
    "input,alpha2",
    [("uk", "UK"), ("el", "EL"), ("de", "DE"), ("Germany", "DE"), ("Greece", "EL")],
)
def test_conversion_to_eu_alpha2(input, alpha2):
    assert convert_country_code(input, "alpha2_eu") == alpha2
