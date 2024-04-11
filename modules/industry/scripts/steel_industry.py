from typing import Optional

import eurocalliopelib.utils as ec_utils
import pandas as pd
from utils import formatting
from utils import jrc_idees_parser as jrc

CAT_NAME_STEEL = "Iron and steel"

H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
HDRI_CONSUMPTION = 0.0116  # H-DRI: 135kWh_e/t = 0.0116ktoe/kt


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
    path_jrc_energy: str,
    path_jrc_production: str,
    path_output: Optional[str] = None,
) -> pd.DataFrame:
    """Execute the data processing pipeline for the "Iron and steel" sub-sector.

    Args:
        year_range (list): years to include and interpolate.
        path_energy_balances (str): country energy balances (usually from eurostat).
        path_cat_names (str): eurostat category mapping file.
        path_carrier_names (str): eurostat carrier name mapping file.
        path_jrc_energy (str): jrc country-specific industrial energy demand file.
        path_jrc_production (str): jrc country-specific industrial production file.
        path_output (str): location of steel demand output file.

    Returns:
        pd.DataFrame: dataframe with steel demand per country.
    """
    # -------------------------------------------------------------------------
    # Prepare data files
    # -------------------------------------------------------------------------
    energy_balances_df = pd.read_csv(
        path_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True
    )
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    energy_df = pd.read_csv(path_jrc_energy, index_col=[0, 1, 2, 3, 4, 5, 6])
    prod_df = pd.read_csv(path_jrc_production, index_col=[0, 1, 2, 3])
    # Ensure dataframes only have data specific to this industry
    cat_names_df = cat_names_df[cat_names_df["jrc_idees"] == CAT_NAME_STEEL]
    energy_df = energy_df.xs(CAT_NAME_STEEL, level="cat_name", drop_level=False)
    prod_df = prod_df.xs(CAT_NAME_STEEL, level="cat_name", drop_level=False)

    # -------------------------------------------------------------------------
    # Process data
    # -------------------------------------------------------------------------
    steel_energy_consumption = process_steel_energy_consumption(
        energy_df, prod_df, cnf_steel
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
    jrc_energy_df: pd.DataFrame, jrc_prod_df: pd.DataFrame, cnf_steel: dict
) -> pd.DataFrame:
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
        jrc_energy_df (pd.DataFrame): jrc country-specific steel energy demand.
        jrc_prod_df (pd.DataFrame): jrc country-specific steel production.
        cnf_steel (dict): configuration for the steel industry.

    Returns:
        pd.DataFrame: processed dataframe with the expected steel energy consumption.
    """

    # sintering/pelletising
    sintering_specific_consumption = jrc.get_specific_electricity_consumption(
        "Integrated steelworks",
        "Steel: Sinter/Pellet making",
        jrc_energy_df,
        jrc_prod_df,
    )
    # smelters
    eaf_smelting_specific_consumption = jrc.get_specific_electricity_consumption(
        "Electric arc", "Steel: Smelters", jrc_energy_df, jrc_prod_df
    )
    # EAF
    eaf_specific_consumption = jrc.get_specific_electricity_consumption(
        "Electric arc", "Steel: Electric arc", jrc_energy_df, jrc_prod_df
    )
    # Rolling & refining
    refining_specific_consumption = jrc.get_specific_electricity_consumption(
        "Electric arc",
        "Steel: Furnaces, Refining and Rolling",
        jrc_energy_df,
        jrc_prod_df,
    )
    finishing_specific_consumption = jrc.get_specific_electricity_consumption(
        "Electric arc", "Steel: Products finishing", jrc_energy_df, jrc_prod_df
    )
    # Auxiliaries (lighting, motors, etc.)
    auxiliary_specific_consumption = jrc.get_auxiliary_electricity_consumption(
        "Electric arc", jrc_energy_df, jrc_prod_df
    )

    # Total electricity consumption
    #   If the country produces steel from Iron ore (assuming 50% recycling):
    #   sintering/pelletizing * iron_ore_% + smelting * recycled_steel_% + H-DRI + EAF + refining/rolling + finishing + auxiliaries
    #   If the country only recycles steel:
    #   smelting + EAF + refining/rolling + finishing + auxiliaries

    recycled_steel_share = cnf_steel["recycled-steel-share"]

    total_specific_consumption = (
        sintering_specific_consumption.mul(1 - recycled_steel_share)
        .add(
            eaf_smelting_specific_consumption
            # if no sintering, this country/year recycles 100% of steel
            .where(sintering_specific_consumption == 0)
            # if there is sintering, update smelting consumption to equal our assumed 2050 recycling rate
            # and add weighted H-DRI consumption to process the remaining iron ore
            .fillna(
                eaf_smelting_specific_consumption.mul(recycled_steel_share).add(
                    HDRI_CONSUMPTION
                )
            )
        )
        .add(eaf_specific_consumption)
        .add(refining_specific_consumption)
        .add(finishing_specific_consumption)
        .add(auxiliary_specific_consumption)
    )
    # In case our model now says a country does produce steel,
    # we give them the average of energy consumption of all other countries
    total_specific_consumption = (
        total_specific_consumption.where(total_specific_consumption > 0)
        .fillna(total_specific_consumption.mean())
        .assign(carrier="electricity")
        .set_index("carrier", append=True)
    )

    # Hydrogen consumption for H-DRI, only for those country/year combinations that handle iron ore
    # and don't recycle all their steel
    h2_specific_consumption = H2_LHV_KTOE * _get_h2_to_steel(recycled_steel_share)
    total_specific_h2_consumption = (
        total_specific_consumption.where(sintering_specific_consumption > 0)
        .fillna(0)
        .where(lambda x: x == 0)
        .fillna(h2_specific_consumption)
        .rename(index={"electricity": "hydrogen"})
    )
    total_specific_consumption = total_specific_consumption.append(
        total_specific_h2_consumption
    )

    # Space heat
    space_heat_specific_demand = (
        jrc_energy_df.xs(("demand", "Electric arc", "Low enthalpy heat"))
        .div(jrc_prod_df.xs("Electric arc").droplevel("unit"))
        .assign(carrier="space_heat")
        .set_index("carrier", append=True)
        .sum(level=total_specific_consumption.index.names)
        .rename(index={"ktoe": "ktoe/kt"})
    )
    total_specific_consumption = total_specific_consumption.append(
        space_heat_specific_demand
    )

    steel_consumption = total_specific_consumption.mul(
        jrc_prod_df.xs("Iron and steel", level="cat_name").sum(level="country_code"),
        level="country_code",
    ).rename(index={"ktoe/kt": "ktoe"})

    return steel_consumption


if __name__ == "__main__":
    get_steel_demand_df(
        year_range=snakemake.params.year_range,
        cnf_steel=snakemake.params.cnf_steel,
        path_energy_balances=snakemake.input.path_energy_balances,
        path_cat_names=snakemake.input.path_cat_names,
        path_carrier_names=snakemake.input.path_carrier_names,
        path_jrc_energy=snakemake.input.path_jrc_energy,
        path_jrc_production=snakemake.input.path_jrc_production,
        path_output=snakemake.output.path_output,
    )
