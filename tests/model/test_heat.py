def test_cop(cop_timeseries):

    # Sanity check that there is a. higher COP in summer than winter, b. no COP < 1 (worse than direct electrical heating)
    cop_monthly = cop_timeseries.groupby(cop_timeseries.index.month).mean()
    cop_winter = cop_monthly.loc[12, 1, 2].mean()
    cop_summer = cop_monthly.loc[6, 7, 8].mean()
    assert (
        cop_summer > cop_winter
    ).all(), "Found higher heating COP values in winter than in summer."
    assert (
        cop_timeseries >= 1
    ).all(), "Found improbably low heat pump COP values (< 1)."


def test_unscaled_heat(unscaled_space_heat_timeseries):

    # Sanity check that there is more space heating demand in winter than summer
    hourly_space_monthly = unscaled_space_heat_timeseries.groupby("time.month").sum()
    hourly_space_winter = hourly_space_monthly.sel(month=[12, 1, 2]).sum("month")
    hourly_space_summer = hourly_space_monthly.sel(month=[6, 7, 8]).sum("month")
    assert (hourly_space_winter > hourly_space_summer).all()
