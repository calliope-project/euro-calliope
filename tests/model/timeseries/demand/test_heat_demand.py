import pytest


@pytest.mark.parametrize("demand", ["heat_demand", "electrified_heat_demand"])
def test_heat_demand_sign(request, demand):
    demand_df = request.getfixturevalue(demand)
    assert (demand_df.stack() >= 0).all(), "Found positive heat demand."


@pytest.mark.parametrize("demand", ["heat_demand", "electrified_heat_demand"])
def test_heat_demand_between_seasons(request, demand, location):
    """Expecting more heating demand in winter than summer."""
    demand_df = request.getfixturevalue(demand)
    demand_monthly = demand_df[location].abs().groupby(demand_df.index.month).sum()
    demand_winter = demand_monthly.loc[[12, 1, 2]].sum()
    demand_summer = demand_monthly.loc[[6, 7, 8]].sum()
    assert (
        demand_winter > demand_summer
    ).all(), "Found higher heat demand in summer than in winter."


def test_electrified_heat_vs_heat_demand(
    heat_demand, electrified_heat_demand, location
):
    """Expecting heat demand to always be the same as or _higher_ than electrified heat demand, due to heating techs having a COP >= 1."""
    assert (
        heat_demand[location] >= electrified_heat_demand[location]
    ).all(), "Found higher electrified heat demand than final heat demand."
