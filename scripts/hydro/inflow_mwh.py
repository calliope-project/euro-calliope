import numpy as np
import pandas as pd
import xarray as xr
from scipy.optimize import minimize


def determine_energy_inflow(path_to_stations_with_water_inflow, path_to_generation, year,
                            max_capacity_factor, path_to_output):
    plants_with_inflow_m3 = xr.open_dataset(path_to_stations_with_water_inflow)
    annual_national_generation_mwh = read_generation(path_to_generation, year)

    inflow_MWh = energy_inflow(plants_with_inflow_m3, annual_national_generation_mwh, max_capacity_factor)
    xr.merge([plants_with_inflow_m3, inflow_MWh]).to_netcdf(path_to_output)


def read_generation(path_to_generation, year):
    return (
        pd
        .read_csv(path_to_generation, index_col=[0, 1])
        .loc[:, "generation_gwh"]
        .rename("generation")
        .xs(year, level="year")
        .mul(1000) # from GWh to MWh
    )


def energy_inflow(plants_with_inflow_m3, annual_national_generation_mwh, max_capacity_factor):
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
    total generation in the year matches the known annual generation. I am still
    capping the time series because of numerical range issues: assume capacity
    factors can be as small as 1e-3, than a capacity factor of 1e1 means a numerical
    range of 10e5 which is large. Uncapped, capacity factors are up to 1e2.
    """
    # ASSUME national fixed capacity / annual generation ratio for run of river and dammed hydro
    # ASSUME dammed hydro never spills water
    # ASSUME dammed hydro inflow is limited
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
            inflow_MWh.loc[{"id": plant_id}] = water_to_capped_energy_inflow(
                inflow_m3=plant.inflow_m3,
                annual_generation=annual_generation_MWh.loc[plant_id],
                cap=plant.installed_capacity_MW.item()
            )
        elif plant.type.item() == "HDAM":
            inflow_MWh.loc[{"id": plant_id}] = water_to_capped_energy_inflow(
                inflow_m3=plant.inflow_m3,
                annual_generation=annual_generation_MWh.loc[plant_id],
                cap=plant.installed_capacity_MW.item() * max_capacity_factor
            )
    return inflow_MWh


def water_to_capped_energy_inflow(inflow_m3, annual_generation, cap):
    assert annual_generation <= cap * 8760

    def generation(scaling_factor):
        generation = inflow_m3 * scaling_factor
        generation[generation > cap] = cap
        return generation

    def residual(scaling_factor):
        return abs(generation(scaling_factor).sum() - annual_generation)

    x0 = annual_generation / inflow_m3.sum()
    res = minimize(residual, x0, method='nelder-mead', options={'xtol': 1e-10})
    assert res.success, print(res)
    assert res.fun < 1, print(res) # error smaller than 1 MWh

    return generation(res.x)


def allocate_generation_to_plant(plants, annual_national_generation_mwh):
    inflows = plants.inflow_m3.to_pandas().fillna(0)
    plants = plants[["country_code", "installed_capacity_MW"]].to_dataframe()

    # Capacity share is scaled to account for erroneous zero timesteps in the inflow data
    inflow_count = inflows.where(inflows > 1).count(axis=1) / inflows.count(axis=1)
    capacity_share = plants.assign(
        scaled_capacity=plants.installed_capacity_MW.multiply(inflow_count)
    ).groupby("country_code")["scaled_capacity"].transform(lambda x: x / x.sum())

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
        max_capacity_factor=snakemake.params.max_capacity_factor,
        path_to_output=snakemake.output[0]
    )
