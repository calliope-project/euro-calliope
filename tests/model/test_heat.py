def test_cop_monthly(cop_timeseries):
    """Expecting higher COP in summer than winter."""
    cop_monthly = cop_timeseries.groupby(cop_timeseries.index.month).mean()
    cop_winter = cop_monthly.loc[[12, 1, 2]].mean()
    cop_summer = cop_monthly.loc[[6, 7, 8]].mean()
    assert (
        cop_summer > cop_winter
    ).all(), "Found higher heating COP values in winter than in summer."


def test_cop_all_at_least_1(cop_timeseries):
    """Expecting no COP < 1 (worse than direct electrical heating)"""
    assert (
        cop_timeseries.stack() >= 1
    ).all(), "Found improbably low heat pump COP values (< 1)."


def test_unscaled_heat(unscaled_space_heat_timeseries):
    """Expecting more space heating demand in winter than summer."""
    # xarray to pandas as we need a newer xarray version (with flox) for `groupby` performance.
    df = unscaled_space_heat_timeseries.space_heat.to_series()
    hourly_space_monthly = df.groupby(df.index.get_level_values("time").month).sum()
    hourly_space_winter = hourly_space_monthly.loc[[12, 1, 2]].sum()
    hourly_space_summer = hourly_space_monthly.loc[[6, 7, 8]].sum()
    assert (hourly_space_winter > hourly_space_summer).all()
