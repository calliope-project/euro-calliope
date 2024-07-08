def test_heat_demand_sign(heat_demand):
    assert (heat_demand.stack() <= 0).all(), "Found positive heat demand."


def test_historic_electrified_heat_sign(historic_electrified_heat):
    assert (
        historic_electrified_heat.stack() >= 0
    ).all(), "Found negative historic electrified heat demand."


def test_heat_demand_between_seasons(heat_demand, location):
    """Expecting more heating demand in winter than summer."""
    demand_monthly = heat_demand[location].abs().groupby(heat_demand.index.month).sum()
    demand_winter = demand_monthly.loc[[12, 1, 2]].sum()
    demand_summer = demand_monthly.loc[[6, 7, 8]].sum()
    assert (
        demand_winter > demand_summer
    ).all(), "Found higher heat demand in summer than in winter."
