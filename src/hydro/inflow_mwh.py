import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import minimize


def determine_energy_inflow(path_to_stations_with_water_inflow, path_to_generation, year, path_to_output):
    plants_with_inflow_m3 = xr.open_dataset(path_to_stations_with_water_inflow)
    annual_national_generation_mwh = read_generation(path_to_generation, year)

    inflow_MWh = energy_inflow(plants_with_inflow_m3, annual_national_generation_mwh)
    xr.merge([plants_with_inflow_m3, inflow_MWh]).to_netcdf(path_to_output)


def read_generation(path_to_generation, year):
    return pd.read_excel(
        path_to_generation,
        sheet_name="Gen 10-18",
        usecols="A,D:G",
        index_col=[0, 1, 2, 3],
        names=["country_code", "product", "year", "type", "generation"]
    )["generation"].to_xarray().sel(type="ONG", year=year, product="Renewable hydropower").to_series() * 1000 # from GWh to MWh


def energy_inflow(plants_with_inflow_m3, annual_national_generation_mwh):
    """Generate hydro power time series based on unscaled water time series.

    As inputs I am using the water inflow time series which is correct in its
    dynamic behavior, but not in its magnitude. To reach time series of usable energy
    with correct magnitude, I am using two different approaches for run of river
    and for dammed hydropower.

    For run of river, I am scaling the water inflow in a way, that the total
    generation in the year matches the known annual generation without ever
    exceeding the installed capacity of the station.

    For dammed hydro electricity, I am assuming the station is never spilling
    any water and thus, I am scaling the water inflow in a way, that the
    total generation in the year matches the known annual generation without
    capping the time series.
    """
    # ASSUME national fixed capacity / annual generation ratio for run of river and dammed hydro
    # ASSUME dammed hydro never spills water
    annual_generation_MWh = allocate_generation_to_plant(plants_with_inflow_m3, annual_national_generation_mwh)
    inflow_MWh = xr.DataArray(
        dims=["id", "time"],
        coords={"id": plants_with_inflow_m3["id"], "time": plants_with_inflow_m3["time"]},
        data=np.ones((len(plants_with_inflow_m3.id), len(plants_with_inflow_m3.time))) * np.nan,
        name="inflow_MWh"
    )
    for plant_id in [id.item() for id in plants_with_inflow_m3.id]:
        plant = plants_with_inflow_m3.sel(id=plant_id)
        if plant.type.item() == "HROR":
            inflow_MWh.loc[{"id": plant_id}] = ror_water_to_energy_inflow(
                plant.inflow_m3,
                annual_generation_MWh.loc[plant_id],
                plant.installed_capacity_MW.item()
            )
        elif plant.type.item() == "HDAM":
            ts_m3 = plant.inflow_m3
            inflow_MWh.loc[{"id": plant_id}] = ts_m3 / ts_m3.sum() * annual_generation_MWh.loc[plant_id]
    return inflow_MWh


def ror_water_to_energy_inflow(inflow_m3, annual_generation, installed_capacity):
    assert annual_generation <= installed_capacity * 8760

    def generation(scaling_factor):
        generation = inflow_m3 * scaling_factor
        generation[generation > installed_capacity] = installed_capacity
        return generation

    def residual(scaling_factor):
        return abs(generation(scaling_factor).sum() - annual_generation)

    x0 = annual_generation / inflow_m3.sum()
    res = minimize(residual, x0, method='nelder-mead', options={'xtol': 1e-10})
    assert res.success, print(res)
    assert res.fun < 1, print(res) # error smaller than 1 MWh

    return generation(res.x)


def allocate_generation_to_plant(plants, annual_national_generation_mwh):
    plants = plants[["country_code", "installed_capacity_MW"]].to_dataframe()
    capacity_share = plants.groupby("country_code")["installed_capacity_MW"].transform(lambda x: x / x.sum())
    national_generation = (
        plants.reset_index()
              .merge(annual_national_generation_mwh, on="country_code", how="left", validate="many_to_one")
              .set_index("id") # without index reset and set, merge removes the index
              .loc[:, "generation"]
    )
    return capacity_share * national_generation


if __name__ == "__main__":
    determine_energy_inflow(
        path_to_stations_with_water_inflow=snakemake.input.stations,
        path_to_generation=snakemake.input.generation,
        year=snakemake.params.year,
        path_to_output=snakemake.output[0]
    )
