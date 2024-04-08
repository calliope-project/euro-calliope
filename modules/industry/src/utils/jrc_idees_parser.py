import numpy as np
import pandas as pd


def get_auxiliary_electricity_consumption(process: str, energy_df: pd.DataFrame, prod_df: pd.DataFrame):
        auxiliaries = ['Lighting', 'Air compressors', 'Motor drives', 'Fans and pumps']
        consumption = (
            energy_df
            .xs(('consumption', process))
            .loc[auxiliaries]
            .sum(level=['country_code', 'unit', 'cat_name'])
        )
        specific_consumption = consumption.div(prod_df.loc[process].droplevel('unit'))
        specific_consumption.index = specific_consumption.index.set_levels(['ktoe/kt'], level='unit')
        return specific_consumption.fillna(0)


def get_specific_electricity_consumption(process, subprocess, energy_df, prod_df):
    consumption = energy_df.xs(('consumption', process, subprocess))
    demand = energy_df.xs(('demand', process, subprocess))
    specific_demand = demand.sum(level=['country_code']).div(prod_df.loc[process].droplevel('unit'))
    efficiency = demand.div(consumption)
    electrical_efficiency = (
        efficiency
        .where(efficiency > 0)
        .xs('Electricity', level='carrier_name')
        .T  # have to reorient the array thanks to a NotImplementedError
        .fillna(efficiency.xs('Electricity', level='carrier_name').mean(axis=1))
        .T  # have to reorient the array thanks to a NotImplementedError
        .fillna(efficiency.xs('Electricity', level='carrier_name').mean())
    )

    specific_consumption = specific_demand.div(electrical_efficiency).rename(index={'ktoe': 'ktoe/kt'})
    assert (
        specific_consumption.fillna(-1).droplevel('unit')
        >= specific_demand.fillna(-1)
    ).all().all()

    return specific_consumption.fillna(0)


def get_carrier_demand(carrier: str, all_demand_df: pd.DataFrame, energy_df: pd.DataFrame) -> pd.DataFrame:
    """
    Get demand for a specific carrier, assuming all end use demand that could consume
    that carrier are completely met by that carrier
    """
    energy = energy_df.xs(carrier, level='carrier_name')
    energy_efficiency = energy.xs('demand').div(energy.xs('consumption'))
    # Fill NaNs (where there is demand, but no consumption in that country)
    # with the average efficiency a. from the country, b. from all countries
    energy_efficiency = energy_efficiency.fillna(energy_efficiency.mean())

    return all_demand_df.reindex(energy_efficiency.index).div(energy_efficiency)
