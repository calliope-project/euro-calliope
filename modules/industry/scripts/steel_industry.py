from typing import Optional

import eurocalliopelib.utils as ec_utils
import pandas as pd
import xarray as xr
from utils import formatting
from utils import jrc_idees_parser as jrc

CAT_NAME_STEEL = "Iron and steel"

H2_LHV_KTOE = 0.0333  # 0.0333 TWh/kt LHV
HDRI_CONSUMPTION = 135e-6  # H-DRI: 135kWh_e/t


def _get_h2_to_steel(recycled_steel_share: float) -> float:
    """Get t_h2/t_steel, usually for H-DRI."""
    # ASSUME: conversion factor of 0.05 t_h2/t_steel.
    return (1 - recycled_steel_share) * 0.05


def get_steel_demand_df(
    year_range: list,
    cnf_steel: dict,
    path_energy_balances: str,
    path_cat_names: str,
    path_carrier_names: str,
    path_jrc_industry_energy: str,
    path_jrc_industry_production: str,
    path_output: Optional[str] = None,
) -> pd.DataFrame:
    """Execute the data processing pipeline for the "Iron and steel" sub-sector.

    Args:
        year_range (list): years to include and interpolate.
        path_energy_balances (str): country energy balances (usually from eurostat).
        path_cat_names (str): eurostat category mapping file.
        path_carrier_names (str): eurostat carrier name mapping file.
        path_jrc_industry_energy (str): jrc country-specific industrial energy demand file.
        path_jrc_industry_production (str): jrc country-specific industrial production file.
        path_output (str): location of steel demand output file.

    Returns:
        pd.DataFrame: dataframe with steel demand per country.
    """
    # -------------------------------------------------------------------------
    # Prepare data files
    # -------------------------------------------------------------------------
    # Read data
    energy_balances_df = pd.read_csv(
        path_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True
    )
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    jrc_energy = xr.open_dataset(path_jrc_industry_energy)
    jrc_prod = xr.open_dataset(path_jrc_industry_production)

    # TODO: fix weird naming convention forced by the JRC module.
    jrc_energy = jrc_energy.rename({"jrc-idees-industry-twh": "value"})
    jrc_prod = jrc_prod.rename({"jrc-idees-industry-twh": "value"})

    # Ensure dataframes only have data specific to this industry
    cat_names_df = cat_names_df[cat_names_df["jrc_idees"] == CAT_NAME_STEEL]
    jrc_energy = jrc_energy.sel(cat_name=CAT_NAME_STEEL)
    jrc_prod = jrc_prod.sel(cat_name=CAT_NAME_STEEL)

    # -------------------------------------------------------------------------
    # Process data
    # -------------------------------------------------------------------------
    steel_energy_consumption = process_steel_energy_consumption(
        jrc_energy, jrc_prod, cnf_steel
    )

    # -------------------------------------------------------------------------
    # Format the final output
    # -------------------------------------------------------------------------
    steel_energy_consumption.columns = steel_energy_consumption.columns.astype(
        int
    ).rename("year")
    filled_consumption_df = formatting.fill_missing_data(
        energy_balances_df,
        cat_names_df,
        carrier_names_df,
        steel_energy_consumption,
        year_range,
    )

    units = filled_consumption_df.index.get_level_values("unit")
    filled_consumption_df.loc[units == "ktoe"] = filled_consumption_df.loc[
        units == "ktoe"
    ].apply(ec_utils.ktoe_to_twh)
    filled_consumption_df = filled_consumption_df.rename({"ktoe": "twh"}, level="unit")
    filled_consumption_df.index = filled_consumption_df.index.set_names(
        "subsector", level="cat_name"
    )
    filled_consumption_df = filled_consumption_df.stack()

    if path_output is not None:
        filled_consumption_df.reorder_levels(formatting.LEVEL_ORDER).to_csv(path_output)

    return filled_consumption_df


def process_steel_energy_consumption(
    jrc_energy: xr.Dataset, jrc_prod: xr.Dataset, cnf_steel: dict
) -> xr.Dataset:
    """Processing of steel energy demand for different carriers.

    Calculates energy consumption in the iron and steel industry based on expected
    change in processes to avoid fossil feedstocks. All process specific energy consumption
    (energy/t_steel) is based on the Electric Arc process (EAF), except sintering, which
    will be required for iron ore processed using H-DRI, but is not required by EAF.

    This function does the following:
    1. Finds all the specific consumption values by getting
        a. process energy demand / produced steel => specific demand
        b. process electrical demand / electrical consumption => electrical efficiency
        c. specific demand / electricial efficiency => specific electricity consumption
    2. Gets total process specific electricity consumption by adding specific consumptions
    for direct electric processes, EAF, H-DRI, smelting, sintering, refining, and finishing
    3. Gets specific hydrogen consumption for all countries that will process iron ore
    4. Gets specific space heat demand based on demand associated with EAF plants
    5. Gets total demand for electricity, hydrogen, and space heat by multiplying specific
    demand by total steel production (by both EAF and BF-BOF routes).

    Args:
        jrc_energy_df (xr.Dataset): jrc country-specific steel energy demand.
        jrc_prod_df (xr.Dataset): jrc country-specific steel production.
        cnf_steel (dict): configuration for the steel industry.

    Returns:
        xr.Dataset: processed dataframe with the expected steel energy consumption.
    """
    # Gather relevant industrial processes
    sintering_intensity = jrc.get_subsection_final_intensity(
        "Integrated steelworks",
        "Steel: Sinter/Pellet making",
        "Integrated steelworks",
        "Electricity",
        jrc_energy,
        jrc_prod,
        fill_empty=True,
    )

    eaf_smelting_intensity = jrc.get_subsection_final_intensity(
        "Electric arc",
        "Steel: Smelters",
        "Electric arc",
        "Electricity",
        jrc_energy,
        jrc_prod,
        fill_empty=True,
    )
    eaf_intensity = jrc.get_subsection_final_intensity(
        "Electric arc",
        "Steel: Electric arc",
        "Electric arc",
        "Electricity",
        jrc_energy,
        jrc_prod,
        fill_empty=True,
    )
    refining_intensity = jrc.get_subsection_final_intensity(
        "Electric arc",
        "Steel: Furnaces, Refining and Rolling",
        "Electric arc",
        "Electricity",
        jrc_energy,
        jrc_prod,
        fill_empty=True,
    )
    finishing_intensity = jrc.get_subsection_final_intensity(
        "Electric arc",
        "Steel: Products finishing",
        "Electric arc",
        "Electricity",
        jrc_energy,
        jrc_prod,
        fill_empty=True,
    )
    auxiliary_intensity = jrc.get_auxiliary_electric_final_intensity(
        "Electric arc", "Electric arc", jrc_energy, jrc_prod, fill_empty=True
    )

    # Total electricity consumption:
    #   If the country produces steel from Iron ore (smelting):
    #   sintering/pelletizing * iron_ore_% + smelting * recycled_steel_% + H-DRI + EAF + refining/rolling + finishing + auxiliaries
    #   If the country only recycles steel:
    #   smelting + EAF + refining/rolling + finishing + auxiliaries
    recycled_share = cnf_steel["recycled-steel-share"]

    # Limit sintering by the share non-recycled steel
    updated_sintering = sintering_intensity * (1 - recycled_share)

    # Update smelting consumption to equal assumed recycling rate
    # and add weighted H-DRI consumption to process the remaining iron ore
    eaf_smelting_recycled = eaf_smelting_intensity * recycled_share + HDRI_CONSUMPTION

    # Countries with no sintering recycle 100% of steel
    updated_eaf_smelting = eaf_smelting_intensity.where(
        sintering_intensity == 0
    ).fillna(eaf_smelting_recycled)

    electric_intensity = (
        updated_sintering
        + updated_eaf_smelting
        + eaf_intensity
        + refining_intensity
        + finishing_intensity
        + auxiliary_intensity
    )
    # In case our model now says a country does produce steel,
    # we give them the average of energy consumption of all other countries per year
    mean_demand_per_year = electric_intensity.mean("country_code")
    electric_intensity = electric_intensity.where(
        electric_intensity > 0, other=mean_demand_per_year
    )
    electric_intensity = electric_intensity.assign_coords(carrier="electricity")

    # Hydrogen consumption for H-DRI:
    # only for country/year that handle iron ore and don't recycle all their steel
    h_dri_h2_intensity = H2_LHV_KTOE * _get_h2_to_steel(recycled_share)

    h2_intensity = electric_intensity.where(sintering_intensity > 0).fillna(0)
    h2_intensity = h2_intensity.where(h2_intensity == 0, h_dri_h2_intensity)
    h2_intensity = h2_intensity.assign_coords(carrier="hydrogen")

    # Space heat
    space_heat_intensity = jrc.get_subsection_useful_intensity(
        "Electric arc", "Low enthalpy heat", "Electric arc", jrc_energy, jrc_prod
    )
    space_heat_intensity = space_heat_intensity.assign_coords(carrier="space_heat")

    total_intensity = xr.concat(
        [electric_intensity, h2_intensity, space_heat_intensity], dim="carrier"
    )

    steel_energy_demand = total_intensity * jrc_prod.sum("produced_material")
    steel_energy_demand["value"].attrs["units"] = "twh/kt"
    breakpoint()
    return steel_energy_demand


if __name__ == "__main__":
    get_steel_demand_df(
        year_range=snakemake.params.year_range,
        cnf_steel=snakemake.params.cnf_steel,
        path_energy_balances=snakemake.input.path_energy_balances,
        path_cat_names=snakemake.input.path_cat_names,
        path_carrier_names=snakemake.input.path_carrier_names,
        path_jrc_industry_energy=snakemake.input.path_jrc_industry_energy,
        path_jrc_industry_production=snakemake.input.path_jrc_industry_production,
        path_output=snakemake.output.path_output,
    )
