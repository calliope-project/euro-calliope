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
    read_load_profiles,
    filter_countries,
    fill_gaps_per_source,
    get_source_choice_per_country
)

THIS_DIR = os.path.dirname(__file__)

class TestLoadHelperFunctions:
    @pytest.mark.parametrize(
        ("min_thresh", "max_tresh", "expected"),
        ([0, 2, [[0, 1], [1, 3]]], [0.5, 2, [[np.nan, 1], [1, 3]]], [0.5, 1.5, [[np.nan, 1], [np.nan, 3]]])
    )
    def test_filter_outliers(self, min_thresh, max_tresh, expected):
        input_load = pd.DataFrame([[0, 1], [1, 3]], columns=["foo", "bar"])
        config = {
            "outlier-data-thresholds": {
                "relative-to-mean-min": min_thresh, "relative-to-mean-max": max_tresh
            }
        }
        filtered = filter_outliers(input_load, config)
        assert np.allclose(expected, filtered.values, equal_nan=True)

    @pytest.mark.parametrize(
        ("interpolate_hours", "expected"),
        ([0, [[np.nan, 0], [0, np.nan], [1, np.nan] , [2, np.nan], [np.nan, 4]]],
         [1, [[0, 0], [0, 1], [1, np.nan] , [2, 3], [2, 4]]],
         [2, [[0, 0], [0, 1], [1, 2] , [2, 3], [2, 4]]])
    )
    def test_interpolate_gaps(self, interpolate_hours, expected):
        nan_data = [[np.nan, 0], [0, np.nan], [1, np.nan] , [2, np.nan], [np.nan, 4]]
        input_load = pd.DataFrame(nan_data, columns=["foo", "bar"])
        interpolated = _interpolate_gaps(input_load, interpolate_hours)
        assert np.allclose(expected, interpolated, equal_nan=True)

    @pytest.mark.parametrize(
        ("input_load", "expected"),
        ([[0, 1, 2, np.nan], [0]],
        [[1, 1, np.nan, 3], [2]],
        [[1, 1, 2.0, 3], []])
    )
    def test_get_index_of_missing_data(self, input_load, expected):
        index = pd.date_range(start="2016-12-29 00:00:00", periods=4, freq="D")
        input_series = pd.Series(input_load, index=index)
        missing_index = _get_index_of_missing_data(input_series, 2016)
        assert (index[expected] == missing_index).all()

    @pytest.mark.parametrize(
        ("this_year", "next_year", "input_timeseries", "expected"),
        ([2016, 2017, ["2016-02-28 23:00:00", "2016-02-29 00:00:00"], ["2016-02-28 23:00:00"]],
        [2016, 2017, ["2016-02-28 22:00:00", "2016-02-28 23:00:00"], ["2016-02-28 22:00:00", "2016-02-28 23:00:00"]],
        [2016, 2020, ["2016-02-28 22:00:00", "2016-02-28 23:00:00"], ["2016-02-28 22:00:00", "2016-02-28 23:00:00"]],
        [2016, 2020, ["2016-02-28 23:00:00", "2016-02-29 00:00:00"], ["2016-02-28 23:00:00", "2016-02-29 00:00:00"]],
        [2015, 2016, ["2015-02-28 22:00:00", "2015-02-28 23:00:00"], ["2015-02-28 22:00:00", "2015-02-28 23:00:00"]])
    )
    def test_ignore_feb_29th(self, this_year, next_year, input_timeseries, expected):
        input_timeseries = pd.to_datetime(input_timeseries)
        expected = pd.to_datetime(expected)
        filtered_timeseries = _ignore_feb_29th(this_year, next_year, input_timeseries)
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
        ("input_load", "expected"),
        ([[[0, 1], [np.nan, 3]], [0]],
        [[[1, np.nan], [np.nan, 3]], [0, 1]],
        [[[0, 1], [0, 3]], [0]])
    )
    def test_countries_with_missing_data_in_model_year(self, input_load, expected):
        filtered = _countries_with_missing_data_in_model_year(pd.DataFrame(input_load))
        assert np.allclose(np.array(expected), filtered, equal_nan=True)

class TestLoadDummyData:
    @staticmethod
    def data_quality_config(
        interpolate_hours=3,
        acceptable_year_diff_for_gap_filling=2,
        fill_29th_feb_from_28th=True,
        data_source_priority_order=["foo", "bar"]
    ):
        return {
            "outlier-data-thresholds": {
                "relative-to-mean-min": 0.25,
                "relative-to-mean-max": 2,
            },
            "interpolate-hours": interpolate_hours,
            "acceptable-year-diff-for-gap-filling": acceptable_year_diff_for_gap_filling,
            "fill-29th-feb-from-28th": fill_29th_feb_from_28th,
            "data-source-priority-order": data_source_priority_order
        }

    def runner(self, data_quality_config):
        path_to_raw_load = os.path.join(THIS_DIR, "..", "resources", "national", "dummy_load.csv")
        data_sources = data_quality_config["data-source-priority-order"]
        countries = ["ALB", "DEU"]
        year = 2016

        raw_load = read_load_profiles(path_to_raw_load, data_sources)
        filtered_load = filter_countries(raw_load, countries)
        filtered_load = filter_outliers(filtered_load, data_quality_config)
        gap_filled_load = pd.concat(
            fill_gaps_per_source(filtered_load, year, data_quality_config, source)
            for source in data_sources
        )
        return get_source_choice_per_country(
            filtered_load.loc[str(year)],
            gap_filled_load,
            data_sources
        )

    def test_success_with_working_data_quality_config(self):
        load = self.runner(self.data_quality_config())
        assert load.isnull().sum().sum() == 0

    @pytest.mark.parametrize(
        ("config_update", "ALB_in_err", "DEU_in_err", "err_date"),
        ([{"interpolate_hours": 1}, True, False, "2016-02-26"],
         [{"acceptable_year_diff_for_gap_filling": 1}, False, True, "2016-02-27"],
         [{"fill_29th_feb_from_28th": False}, True, False, "2016-02-29"],
         [{"data_source_priority_order": ["foo"]}, False, True, "2016-02-27"])
    )
    def test_fail_on_too_low_interpolate(self, config_update, ALB_in_err, DEU_in_err, err_date):
        config = self.data_quality_config(**config_update)
        with pytest.raises(AssertionError) as excinfo:
            self.runner(config)
        for country, country_error in {"ALB": ALB_in_err, "DEU": DEU_in_err}.items():
            if country_error is True:
                assert country in str(excinfo.value)
            else:
                assert country not in str(excinfo.value)
        assert err_date in str(excinfo.value)
