"""Generate Calliope load time series."""

import pandas as pd
import xarray as xr
import geopandas as gpd
import math


def load(path_units_with_population_share: str,
         path_units_share_of_nat_demand_per_sector: str,
         path_nat_ind_elec_demand: str,
         path_nat_elec_load: str,
         path_result: str,
         scaling_factor_power: float) -> None:
    """
    Generate load time series for every location.

    Parameters:
        path_units_with_population_share (str):
            Units with population share. Used to split/distribute the national residential elec load to units.
        path_units_share_of_nat_demand_per_sector (str):
            Units' fractions of national electricity demand within industrial sectors, computed from emissions.
        path_nat_ind_elec_demand (str):
            Nations' electricity demand for each industrial sector.
        path_nat_elec_load (str):
            Nations' total electricity load (industrial and residential).
        path_result (str):
            Output filepath for .csv file with total electricity load for specified resolution.
        scaling_factor_power (float):
            Scaling for Calliope's solving algorithm.
    """

    units = gpd.read_file(path_units_with_population_share).set_index("id")
    units_share_of_nat_demand_per_sector = pd.read_csv(path_units_share_of_nat_demand_per_sector, index_col=0)
    nat_ind_elec_demand = pd.read_csv(path_nat_ind_elec_demand, index_col=0)
    nat_elec_load = pd.read_csv(path_nat_elec_load, index_col=0, parse_dates=True)

    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        nat_ind_elec_demand = pd.DataFrame(nat_ind_elec_demand.sum(axis=1).rename("EUR"))
        nat_elec_load = pd.DataFrame(nat_elec_load.sum(axis=1).rename("EUR"))

    ind_demand_per_sector_and_unit = compute_ind_demand_per_sector_and_unit(units_share_of_nat_demand_per_sector,
                                                                            nat_ind_elec_demand,
                                                                            units)
    ind_load_per_unit = compute_ind_load_per_unit(ind_demand_per_sector_and_unit,
                                                  nat_elec_load)
    res_load_per_unit = compute_res_load_per_unit(ind_load_per_unit,
                                                  nat_elec_load,
                                                  units)

    load_per_unit = (ind_load_per_unit + res_load_per_unit) * (-1) * scaling_factor_power
    load_per_unit = load_per_unit.to_dataframe("load_per_unit")["load_per_unit"].unstack()
    load_per_unit.columns = load_per_unit.columns.str.replace(".", "-", regex=False)

    assert math.isclose(
        load_per_unit.sum().sum() * (-1) / scaling_factor_power,
        nat_elec_load.reindex(columns=units.country_code.unique()).sum().sum()
    )

    load_per_unit.tz_convert(None).to_csv(path_result)


def compute_ind_demand_per_sector_and_unit(units_share_of_nat_demand_per_sector: pd.DataFrame,
                                           nat_ind_elec_demand: pd.DataFrame,
                                           units: gpd.GeoDataFrame) -> xr.DataArray:

    """
    Computes industrial electricity demand of each industrial sector in each unit by multiplying the sectors' national
    demands with the unit's share of the national emissions in that sector.

    Parameters:
        units (gpd.GeoDataFrame):
            Units with population share. Used for list of country codes.
        units_share_of_nat_demand_per_sector (pd.DataFrame):
            Units' fractions of national electricity demand within industrial sectors, computed from emissions.
        nat_ind_elec_demand (pd.DataFrame):
            Nations' electricity demand for each industrial sector.
    Returns:
        ind_demand_per_sector_and_unit (pd.DataFrame):
            Industrial demand for all units and all industrial sectors.
    """

    units_share_of_nat_demand_per_sector = (
        units_share_of_nat_demand_per_sector
        .stack()
        .rename_axis(index=["industry_sector", "units"])
        .to_xarray()
        .assign_coords(country_code=units["country_code"].rename_axis(index="units"))
    )

    # Inflate nat_ind_elec_demand to have as many columns for each nation as there are units within that nation.
    nat_ind_elec_demand = (
        nat_ind_elec_demand
        .reindex(columns=units["country_code"].values)
    )
    nat_ind_elec_demand.columns = units.index

    nat_ind_elec_demand = (
        nat_ind_elec_demand
        .stack()
        .rename_axis(index=["industry_sector", "units"])
        .to_xarray()
    )

    # ASSUME that national demand can be distributed to units proportional to their emissions
    # Note that nat_ind_elec_demand includes all industrial sectors as defined in the EUROSTAT balances (not transport
    # and other). units_share_of_nat_demand_per_sector includes all industrial sectors (inlcuding transport and other)
    # that have at least one installation in any country in the scope. ind_demand_per_sector_and_unit only has those
    # industrial sectors, that occur in both dataframes (inner joint/multiplication).
    ind_demand_per_sector_and_unit = nat_ind_elec_demand * units_share_of_nat_demand_per_sector

    return ind_demand_per_sector_and_unit


def compute_ind_load_per_unit(ind_demand_per_sector_and_unit: xr.DataArray,
                              nat_elec_load: pd.DataFrame) -> xr.DataArray:
    """
    Computes industrial electricity load by uniformly distributing the yearly demand over all timesteps. First find
    loadcurve for each unit and industrial sector, then sum over industrial sectors to obtain total industrial load for
    each unit.

    Parameters:
        ind_demand_per_sector_and_unit (pd.DataFrame):
            Industrial demand for all units and all industrial sectors.
        nat_elec_load (pd.DataFrame):
            Nations' total electricity load (industrial and residential).
    Returns:
        ind_load_per_unit (pd.DataFrame):
            Total industrial load curves for each unit.

    """

    timesteps = nat_elec_load.index
    # ASSUME flat industry load profiles
    # Includes all industrial sectors (not transport and others) that have at least one installation in any country of
    # the scope. If returning ind_load_per_sector_and_unit and working with the sector's load timeseries, consider
    # including 0-timeseries for the remaining sectors.
    ind_load_per_sector_and_unit = ind_demand_per_sector_and_unit.expand_dims(utc_timestamp=timesteps) / len(timesteps)
    ind_load_per_unit = ind_load_per_sector_and_unit.sum("industry_sector")

    return ind_load_per_unit


def compute_res_load_per_unit(ind_load_per_unit: xr.DataArray,
                              nat_elec_load: pd.DataFrame,
                              units: gpd.GeoDataFrame) -> xr.DataArray:
    """
    Computes residential electricity load for each unit. First, compute national residential load by substracting
    industrial load from national load. For this, sum industrial load over all units within a nation. Obtain
    residential load of each unit by splitting national residential load according to population share of a unit
    within its nation.

    Parameters:
        ind_demand_per_sector_and_unit (pd.DataFrame):
            Industrial demand for all units and all industrial sectors.
        nat_elec_load (pd.DataFrame):
            Nations' total electricity load (industrial and residential).
    Returns:
        ind_load_per_unit (pd.DataFrame):
            Total industrial load curves for each unit.
    """

    nat_elec_load = nat_elec_load.to_xarray().to_array("country_code")
    nat_res_elec_load = nat_elec_load - ind_load_per_unit.groupby("country_code").sum("units")

    nat_res_elec_load = (
        nat_res_elec_load
        .to_dataframe("nat_res_elec_load")["nat_res_elec_load"]
        .unstack(level=0)
        .reindex(columns=units["country_code"].values))
    nat_res_elec_load.columns = units.index

    res_load_per_unit = nat_res_elec_load * units["population_share"]

    return res_load_per_unit.to_xarray().to_array("units")


if __name__ == '__main__':
    load(
        path_units_with_population_share=snakemake.input.units_with_population_share,
        path_units_share_of_nat_demand_per_sector=snakemake.input.units_share_of_nat_demand_per_sector,
        path_nat_ind_elec_demand=snakemake.input.nat_ind_elec_demand,
        path_nat_elec_load=snakemake.input.national_load,
        path_result=snakemake.output[0],
        scaling_factor_power=snakemake.params.scaling_factor_power
    )
