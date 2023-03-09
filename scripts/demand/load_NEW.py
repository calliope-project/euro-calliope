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
    """Generate load time series for every location."""

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

    assert math.isclose(
        load_per_unit.sum().sum() * (-1) / scaling_factor_power,
        nat_elec_load.reindex(columns=units.country_code.unique()).sum().sum()
    )

    load_per_unit.tz_convert(None).to_csv(path_result)


def compute_ind_demand_per_sector_and_unit(units_share_of_nat_demand_per_sector: pd.DataFrame,
                                           nat_ind_elec_demand: pd.DataFrame,
                                           units: gpd.GeoDataFrame) -> xr.DataArray:

    units_share_of_nat_demand_per_sector = (
        units_share_of_nat_demand_per_sector
        .stack()
        .rename_axis(index=["industry_sector", "units"])
        .to_xarray()
        .assign_coords(country_code=units["country_code"].rename_axis(index="units"))
    )

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

    ind_demand_per_sector_and_unit = nat_ind_elec_demand * units_share_of_nat_demand_per_sector

    return ind_demand_per_sector_and_unit


def compute_ind_load_per_unit(ind_demand_per_sector_and_unit: xr.DataArray,
                              nat_elec_load: pd.DataFrame) -> xr.DataArray:

    timesteps = nat_elec_load.index
    # ASSUME flat industry load profiles
    ind_load_per_sector_and_unit = ind_demand_per_sector_and_unit.expand_dims(utc_timestamp=timesteps) / len(timesteps)
    ind_load_per_unit = ind_load_per_sector_and_unit.sum("industry_sector")

    return ind_load_per_unit


def compute_res_load_per_unit(ind_load_per_unit: xr.DataArray,
                              nat_elec_load: pd.DataFrame,
                              units: gpd.GeoDataFrame) -> xr.DataArray:

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
