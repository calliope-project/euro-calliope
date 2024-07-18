import pandas as pd
import xarray as xr
from utils import filling
from utils import jrc_idees_parser as jrc

CAT_NAME_CHEMICALS = "Chemicals Industry"


def get_chemicals_demand_df(
    config: dict,
    energy_balances: str,
    cat_names: str,
    carrier_names: str,
    jrc_industry_energy: str,
    jrc_industry_production: str,
    output_path: str,
):
    """Execute the data processing pipeline for the "Chemicals Industry" sub-sector.

    Args:
        config (dict): chemicals industry sector configuration.
        energy_balances (str): country energy balances (usually from eurostat).
        cat_names (str): eurostat category mapping file.
        carrier_names (str): eurostat carrier name mapping file.
        jrc_industry_energy (str): jrc country-specific industrial energy demand file.
        jrc_industry_production (str): jrc country-specific industrial production file.
        output_path (str): location of chemicals demand output file. Defaults to None.
    """
    # Load data
    energy_balances_df = pd.read_csv(
        energy_balances, index_col=[0, 1, 2, 3, 4]
    ).squeeze("columns")
    cat_names_df = pd.read_csv(cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(carrier_names, header=0, index_col=0)
    jrc_energy = xr.open_dataset(jrc_industry_energy)
    jrc_prod = xr.open_dataarray(jrc_industry_production)
    jrc.check_units(jrc_energy, jrc_prod)

    # Ensure dataframes only have data specific to this industry
    cat_names_df = cat_names_df[cat_names_df["jrc_idees"] == CAT_NAME_CHEMICALS]
    jrc_energy = jrc_energy.sel(cat_name=CAT_NAME_CHEMICALS)
    jrc_prod = jrc_prod.sel(cat_name=CAT_NAME_CHEMICALS)

    # Process data
    new_chemicals_demand = transform_jrc_chemicals_subsector_demand(
        jrc_energy, jrc_prod, config
    )
    new_chemicals_demand = filling.fill_missing_countries_years(
        energy_balances_df, cat_names_df, carrier_names_df, new_chemicals_demand
    )
    new_chemicals_demand.to_netcdf(output_path)


def transform_jrc_chemicals_subsector_demand(
    jrc_energy: xr.Dataset, jrc_prod: xr.Dataset, config_steel: dict
) -> xr.Dataset:
    """Processing of chemicals industry energy demand for different carriers.

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
        config_steel (dict): configuration for the steel industry.

    Returns:
        xr.Dataset: processed dataframe with the expected steel energy consumption.
    """
