def test_cop_monthly(cop_timeseries, location):
    """Expecting higher COP in summer than winter."""
    cop_monthly = cop_timeseries[location].groupby(cop_timeseries.index.month).mean()
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


def test_historic_electrified_heat_sign(historic_electrified_heat):
    assert (
        historic_electrified_heat.stack() >= 0
    ).all(), "Found negative historic electrified heat demand."
