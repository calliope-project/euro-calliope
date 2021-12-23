"""Generate Calliope load time series."""
import math

import pandas as pd
import geopandas as gpd


def load(path_to_units, path_to_demand_per_unit, path_to_electricity_load, scaling_factor, path_to_result):
    """Generate load time series for every location."""
    units = gpd.read_file(path_to_units).set_index("id")
    demand_per_unit = pd.read_csv(path_to_demand_per_unit, index_col=0)
    national_load = pd.read_csv(path_to_electricity_load, index_col=0, parse_dates=True)

    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        national_load = pd.DataFrame(national_load.sum(axis=1).rename("EUR"))
    # demand_per_unit.demand_twh_per_year is not necessarily the demand in the year
    # used here and thus its absolute value must be ignored.
    units["industrial_demand"] = demand_per_unit.demand_twh_per_year * demand_per_unit.industrial_demand_fraction
    units["residential_demand"] = demand_per_unit.demand_twh_per_year - units.industrial_demand
    assert not units["industrial_demand"].isna().any()
    units["fraction_of_national_industrial_load"] = units.groupby("country_code").industrial_demand.transform(
        lambda x: x / x.sum()
    ).fillna(0) # if national demand is 0, division by zero
    units["fraction_of_national_residential_load"] = units.groupby("country_code").residential_demand.transform(
        lambda x: x / x.sum()
    )

    national_industrial_load, national_residential_load = split_national_load(national_load, units)
    load_ts = pd.concat(
        [unit_time_series(unit, unit_name, national_industrial_load, national_residential_load, scaling_factor)
         for unit_name, unit in units.iterrows()],
        axis=1
    )
    assert math.isclose(
        load_ts.sum().sum() * (-1) / scaling_factor,
        national_load.reindex(columns=units.country_code.unique()).sum().sum()
    )
    load_ts.tz_convert(None).to_csv(path_to_result)


def split_national_load(national_load, units):
    national_industrial_demand = units.groupby("country_code").industrial_demand.sum() * 1e6 # from TWh to MWh
    industrial_load = pd.DataFrame( # ASSUME flat industry load profiles
        index=national_load.index,
        data=national_industrial_demand.to_dict()
    ).div(len(national_load.index)).reindex(columns=national_load.columns, fill_value=0)
    residential_load = national_load - industrial_load
    return industrial_load, residential_load


def unit_time_series(unit, unit_name, national_industrial_load, national_residential_load, scaling_factor):
    country_code = unit.country_code
    multiplier = unit.fraction_of_national_industrial_load
    unit_industrial_ts = national_industrial_load.loc[:, country_code].copy() * multiplier * (-1) * scaling_factor
    multiplier = unit.fraction_of_national_residential_load
    unit_residential_ts = national_residential_load.loc[:, country_code].copy() * multiplier * (-1) * scaling_factor
    unit_ts = unit_industrial_ts + unit_residential_ts
    unit_ts.name = unit_name.replace(".", "-")
    return unit_ts


if __name__ == '__main__':
    load(
        path_to_units=snakemake.input.units,
        path_to_demand_per_unit=snakemake.input.demand_per_unit,
        path_to_electricity_load=snakemake.input.national_load,
        path_to_result=snakemake.output[0],
        scaling_factor=snakemake.params.scaling_factor
    )
