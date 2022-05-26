"""
Use Eurostat/JRC-IDEES data to scale hourly heat profiles
"""

import pandas as pd
import xarray as xr


def scale_heat_profiles(
    paths_to_heat_profiles,
    path_to_dwelling_ratios,
    path_to_annual_heat_demand,
    resolution,
    reference_year,
    scaling_factor_power,
    path_to_output
):
    all_heat_profiles = xr.concat(
        [xr.open_dataarray(file) for file in paths_to_heat_profiles], dim="time"
    ).sortby("time")

    normalised_heat_profiles = all_heat_profiles / all_heat_profiles.sel(time=str(reference_year)).sum("time")
    annual_heat_demand = xr.open_dataarray(path_to_annual_heat_demand).sel(year=reference_year)

    if resolution == "continental":
        annual_heat_demand = annual_heat_demand.sum("id")

    commercial_profiles = scale_commercial_profiles(normalised_heat_profiles, annual_heat_demand)
    household_profiles = scale_household_profiles(normalised_heat_profiles, annual_heat_demand, path_to_dwelling_ratios)

    (
        (commercial_profiles + household_profiles)
        .reset_coords(drop=True)
        .to_series()
        .unstack("id")
        .mul(-1 * scaling_factor_power)  # MW -> model scale & make negative to act as a calliope input
        .to_csv(path_to_output)
    )


def scale_commercial_profiles(normalised_heat_profiles, annual_heat_demand):
    # Commercial demand is a single profile and a single annual value
    annual_commercial_demand = annual_heat_demand.sel(cat_name="commercial")
    return normalised_heat_profiles.sel(building="COM") * annual_commercial_demand


def scale_household_profiles(normalised_heat_profiles, annual_heat_demand, path_to_dwelling_ratios):
    # Household demand has two profiles, which are merged based on the ratio of single-
    # to multi-family homes before scaling to the total annual household demand
    dwellings_ratios = pd.read_csv(path_to_dwelling_ratios, index_col=0).squeeze()  # MFH/SFH
    annual_household_demand = annual_heat_demand.sel(cat_name="household")
    normalised_heat_profiles = (
        normalised_heat_profiles.sel(building="MFH") * dwellings_ratios.to_xarray() +
        normalised_heat_profiles.sel(building="SFH") * (1 - dwellings_ratios).to_xarray()
    )
    return normalised_heat_profiles * annual_household_demand


if __name__ == "__main__":
    scale_heat_profiles(
        paths_to_heat_profiles=snakemake.input.heat_profiles,
        path_to_dwelling_ratios=snakemake.input.dwelling_ratios,
        path_to_annual_heat_demand=snakemake.input.annual_heat_demand,
        reference_year=snakemake.params.reference_year,
        scaling_factor_power=snakemake.params.scaling_factor_power,
        resolution=snakemake.wildcards.resolution,
        path_to_output=snakemake.output[0]
    )
