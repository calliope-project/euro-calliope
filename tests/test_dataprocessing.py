import pandas as pd


def test_annual_energy_balance(annual_energy_balance):
    assert isinstance(annual_energy_balance, pd.DataFrame)

    # assert annual_energy_balance.columns == [...]

    # assert sum(annual_energy_balance["FOOBAR"]) == BAZ
