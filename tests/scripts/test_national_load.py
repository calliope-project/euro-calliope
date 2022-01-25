import os
import pytest

import pandas as pd
import numpy as np

from scripts.national_load import (
    filter_outliers,
    _interpolate_gaps,
    _fill_29th_feb,
    _countries_with_missing_data_in_model_year,
    _get_index_of_missing_data,
    _ignore_feb_29th,
    clean_load_data
)

THIS_DIR = os.path.dirname(__file__)

class TestLoadHelperFunctions:
    @pytest.fixture
    def foobar_df(self):
        def _foobar_df(foo, bar):
            return pd.DataFrame({"foo": foo, "bar": bar})
        return _foobar_df


    @pytest.mark.parametrize(
        ("min_thresh", "max_tresh", "expected_foo", "expected_bar"), (
            [0, 2, [0, 1], [1, 3]],
            [0.5, 2, [np.nan, 1], [1, 3]],
            [0.75, 1.5, [np.nan, np.nan], [np.nan, 3]]
        )
    )
    def test_filter_outliers(self, min_thresh, max_tresh, expected_foo, expected_bar, foobar_df):
        input_load = foobar_df(
            foo=[0, 1], # mean: 0.5
            bar=[1, 3]  # mean: 2
        )
        config = {
            "outlier-data-thresholds": {
                "relative-to-mean-min": min_thresh, "relative-to-mean-max": max_tresh
            }
        }
        filtered = filter_outliers(input_load, config)
        expected = foobar_df(expected_foo, expected_bar)
        assert np.allclose(expected, filtered.values, equal_nan=True)

    @pytest.mark.parametrize(
        ("interpolate_timesteps", "expected_foo", "expected_bar"),
        ([0, [np.nan, 0, 1, 2, np.nan], [0, np.nan, np.nan, np.nan, 4]],
         [1, [0, 0, 1, 2, 2], [0, 1, np.nan, 3, 4]],
         [2, [0, 0, 1, 2, 2], [0, 1, 2, 3, 4]])
    )
    def test_interpolate_gaps(self, interpolate_timesteps, expected_foo, expected_bar, foobar_df):
        input_load = foobar_df(
            foo=[np.nan, 0, 1, 2, np.nan],
            bar=[0, np.nan, np.nan, np.nan, 4]
        )
        interpolated = _interpolate_gaps(input_load, interpolate_timesteps)
        expected = foobar_df(expected_foo, expected_bar)
        assert np.allclose(expected, interpolated, equal_nan=True)

    @pytest.mark.parametrize(
        ("input_load", "year", "expected"),
        ([[0, 1, 2, np.nan], 2016, [0]],
         [[1, 1, np.nan, 3], 2016, [2]],
         [[1, 1, 2.0, 3], 2016, []],
         [[0, 1, 2, np.nan], 2017, [3]],
         [[1, 1, np.nan, 0], 2017, [3]],
         [[1, 1, 2.0, 3], 2017, []])
    )
    def test_get_index_of_missing_data(self, input_load, year, expected):
        # index = 2016-12-29 to 2017-01-01;
        # `_get_index_of_missing_data` will only return index values for the year of interest
        index = pd.date_range(start="2016-12-29 00:00:00", periods=4, freq="D")
        input_series = pd.Series(input_load, index=index)
        missing_index = _get_index_of_missing_data(input_series, year)
        assert (index[expected] == missing_index).all()

    @pytest.mark.parametrize(
        ("model_year", "next_available_year_of_data", "input_timeseries", "expected"),
        ([2016, 2017, ["2016-02-28 23:00:00", "2016-02-29 00:00:00"], ["2016-02-28 23:00:00"]],
        [2016, 2017, ["2016-02-28 22:00:00", "2016-02-28 23:00:00"], ["2016-02-28 22:00:00", "2016-02-28 23:00:00"]],
        [2016, 2020, ["2016-02-28 22:00:00", "2016-02-28 23:00:00"], ["2016-02-28 22:00:00", "2016-02-28 23:00:00"]],
        [2016, 2020, ["2016-02-28 23:00:00", "2016-02-29 00:00:00"], ["2016-02-28 23:00:00", "2016-02-29 00:00:00"]],
        [2015, 2016, ["2015-02-28 22:00:00", "2015-02-28 23:00:00"], ["2015-02-28 22:00:00", "2015-02-28 23:00:00"]])
    )
    def test_ignore_feb_29th(self, model_year, next_available_year_of_data, input_timeseries, expected):
        input_timeseries = pd.to_datetime(input_timeseries)
        expected = pd.to_datetime(expected)
        filtered_timeseries = _ignore_feb_29th(
            model_year, next_available_year_of_data, input_timeseries
        )
        assert (filtered_timeseries == expected).all()

    @pytest.mark.parametrize(
        ("year", "input_load", "expected"),
        ([2016, {"2016-02-28 23:00:00": 1}, {"2016-02-28 23:00:00": 1}],
         [2017, {"2017-02-28 23:00:00": 1}, {"2017-02-28 23:00:00": 1}],
        [2016,
        {"2016-02-28 22:00:00": 1, "2016-02-28 23:00:00": 2, "2016-02-29 22:00:00": np.nan, "2016-02-29 23:00:00": 3},
        {"2016-02-28 22:00:00": 1, "2016-02-28 23:00:00": 2, "2016-02-29 22:00:00": 1, "2016-02-29 23:00:00": 3}])
    )
    def test_fill_29th_feb(self, year, input_load, expected):
        input_load = pd.Series(input_load)
        input_load.index = pd.to_datetime(input_load.index)
        expected = pd.Series(expected)
        expected.index = pd.to_datetime(expected.index)
        filled_load = _fill_29th_feb(input_load, year)
        assert expected.astype(float).equals(filled_load.astype(float))

    @pytest.mark.parametrize(
        ("input_foo", "input_bar", "expected"),
        ([[0, np.nan], [1, 3], ["foo"]],
         [[1, np.nan], [np.nan, 3], ["foo", "bar"]],
        [[0, 0], [1, 3], ["foo"]])
    )
    def test_countries_with_missing_data_in_model_year(self, input_foo, input_bar, expected, foobar_df):
        input_load = foobar_df(input_foo, input_bar)
        filtered = _countries_with_missing_data_in_model_year(input_load)
        assert filtered.difference(expected).empty

class TestLoadDummyData:
    @pytest.fixture
    def load(self):
        def _load(
            interpolate_timesteps=3,
            acceptable_year_diff_for_gap_filling=2,
            fill_29th_feb_from_28th=True,
            data_source_priority_order=["foo", "bar"]
        ):
            data_quality_config = {
                "outlier-data-thresholds": {
                    "relative-to-mean-min": 0.25,
                    "relative-to-mean-max": 2,
                },
                "max-interpolate-timesteps": interpolate_timesteps,
                "acceptable-year-diff-for-gap-filling": acceptable_year_diff_for_gap_filling,
                "fill-29th-feb-from-28th": fill_29th_feb_from_28th,
                "data-source-priority-order": data_source_priority_order
            }
            path_to_raw_load = os.path.join(THIS_DIR, "..", "resources", "national", "dummy_load.csv")
            countries = ["ALB", "DEU"]
            year = 2016

            return clean_load_data(path_to_raw_load, year, year, data_quality_config, countries)
        return _load

    def test_success_with_working_data_quality_config(self, load):
        assert load().isnull().sum().sum() == 0

    # regex checks if the country names are in the error message, `?=` = `is in`, `?!=` = `not in`
    @pytest.mark.parametrize(
        ("config_update", "countries_in_err_msg", "err_date"),
        ([{"interpolate_timesteps": 1}, r"(?=.*ALB)(?!=.*DEU)", "2016-02-26"],
         [{"acceptable_year_diff_for_gap_filling": 1}, r"(?!=.*ALB)(?=.*DEU)", "2016-02-27"],
         [{"fill_29th_feb_from_28th": False}, r"(?=.*ALB)(?!=.*DEU)", "2016-02-29"],
         [{"data_source_priority_order": ["foo"]}, r"(?!=.*ALB)(?=.*DEU)", "2016-02-27"])
    )
    def test_fail_on_too_low_interpolate(
        self, load, config_update, countries_in_err_msg, err_date
    ):
        with pytest.raises(AssertionError, match=countries_in_err_msg) as excinfo:
            load(**config_update)
        assert err_date in str(excinfo.value)
