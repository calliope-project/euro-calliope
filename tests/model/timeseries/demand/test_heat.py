def test_unscaled_heat(unscaled_space_heat_timeseries):
    """Expecting more space heating demand in winter than summer."""
    # xarray to pandas as we need a newer xarray version (with flox) for `groupby` performance.
    df = unscaled_space_heat_timeseries.space_heat.to_series()
    hourly_space_monthly = df.groupby(df.index.get_level_values("time").month).sum()
    hourly_space_winter = hourly_space_monthly.loc[[12, 1, 2]].sum()
    hourly_space_summer = hourly_space_monthly.loc[[6, 7, 8]].sum()
    assert (hourly_space_winter > hourly_space_summer).all()
