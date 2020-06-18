from datetime import datetime

import pytest


def test_capacity_factors_are_trimmed(capacity_factor_timeseries, config):
    if not config["capacity-factors"]["trim-ninja-timeseries"]:
        pytest.skip("Trim not demanded.")
    year = config["year"]
    first_index = capacity_factor_timeseries.index[0]
    last_index = capacity_factor_timeseries.index[-1]
    assert first_index >= datetime(year=year, month=1, day=1)
    assert last_index < datetime(year=year + 1, month=1, day=1)


def test_capacity_factors_are_floored(capacity_factor_timeseries, config):
    floor = config["capacity-factors"]["min"]
    assert capacity_factor_timeseries[capacity_factor_timeseries > 0].min().min() >= floor


def test_capacity_factors_are_capped(capacity_factor_timeseries, config):
    cap = config["capacity-factors"]["max"]
    assert capacity_factor_timeseries[capacity_factor_timeseries > 0].max().max() <= cap


def test_open_field_vs_rooftop(open_field_pv_capacity_factor_timeseries,
                               rooftop_pv_capacity_factor_timeseries, location):
    of = open_field_pv_capacity_factor_timeseries.loc[:, location].mean()
    rt = rooftop_pv_capacity_factor_timeseries.loc[:, location].mean()
    assert of > rt
