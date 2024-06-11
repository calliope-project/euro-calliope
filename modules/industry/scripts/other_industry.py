from typing import Optional

import pandas as pd
import xarray as xr
from utils import filling
from utils import jrc_idees_parser as jrc


def get_generic_demand(
    non_generic_categories: list,
    generic_config: dict,
    path_energy_balances: str,
    path_cat_names: str,
    path_carrier_names: str,
    path_jrc_industry_energy: str,
    path_jrc_industry_production: str,
    path_output: Optional[str] = None,
) -> xr.DataArray:
    """Execute a generic data processing pipeline for industry categories.

    Args:
        non_generic_categories (list): categories with separate processing (will be ignored).
        generic_config (dict): configuration for generic category processing.
        path_energy_balances (str): country energy balances (usually from eurostat).
        path_cat_names (str): eurostat category mapping file.
        path_carrier_names (str): eurostat carrier name mapping file.
        path_jrc_industry_energy (str): jrc country-specific industrial energy demand file.
        path_jrc_industry_production (str): jrc country-specific industrial production file.
        path_output (str): location of steel demand output file.

    Returns:
        pd.DataFrame: dataframe with industrial demand per country.
    """
    # Load data
    energy_balances_df = pd.read_csv(
        path_energy_balances, index_col=[0, 1, 2, 3, 4]
    ).squeeze("columns")
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    jrc_energy = xr.open_dataset(path_jrc_industry_energy)
    jrc_prod = xr.open_dataarray(path_jrc_industry_production)

    # Remove data from all specifically processed industries
    cat_names_df = cat_names_df[~cat_names_df["jrc_idees"].isin(non_generic_categories)]
    jrc_energy = jrc_energy.drop_sel(cat_name=non_generic_categories)
    jrc_prod = jrc_prod.drop_sel(cat_name=non_generic_categories)

    # Process data:
    # Extract useful dem. -> remove useful dem. from rest -> extract final dem.
    selected_useful = generic_config["useful-demands"]
    other_useful_demand = jrc.convert_subsec_demand_to_carrier(
        jrc_energy, selected_useful
    )

    final_method = generic_config["final-energy-method"]
    jrc_energy = jrc_energy.drop_sel(subsection=selected_useful)

    match final_method:
        case "by priority":
            other_final_demand = transform_final_demand_by_priority(
                jrc_energy, generic_config["final-energy-carriers"]
            )
        case "keep everything":
            other_final_demand = jrc_energy["final"].sum(["section", "subsection"])
            other_final_demand = jrc.standardize(other_final_demand, "twh")
        case _:
            raise ValueError(f"Unsupported final energy method: {final_method}.")

    # Combine and fill missing countries
    other_demand = xr.concat(
        [other_useful_demand, other_final_demand], dim="carrier_name"
    )

    other_demand = filling.fill_missing_countries_years(
        energy_balances_df, cat_names_df, carrier_names_df, other_demand
    )

    other_demand = jrc.standardize(other_demand, "twh")

    if path_output:
        other_demand.to_netcdf(path_output)

    return other_demand


def transform_final_demand_by_priority(
    jrc_energy: xr.Dataset, carrier_priority: list[str]
) -> xr.DataArray:
    """Transform final demand of generic categories by giving priority to certain carriers.

    Steps:
    1. Assume that all demand that could consume a carrier will be met by said carrier.
    2. Drop overlapping consumption so that demand is met by carriers with the given priority.
    3. Combine.

    E.g., if carrier priority is [Electricity, Natural gas, Diesel] then:
    - Electricity: if met exclusively or otherwise, it's final electrical demand.
    - Natural gas: if met exclusively or otherwise, EXCEPT for overlapping cases with Electricity.
    - Diesel: if met exclusively or otherwise, EXCEPT for overlapping cases with all the above.

    Args:
        jrc_energy (xr.Dataset): JRC energy dataset.
        carrier_priority (list[str]): carriers to take in order of priority.

    Returns:
        xr.DataArray: dataset filled with demands for the given carriers.
    """
    carrier_final_dem = {}

    for carrier in carrier_priority:
        dem_replaced = jrc.replace_final_demand_by_carrier(carrier, jrc_energy)
        dem_replaced = dem_replaced.to_dataframe().dropna()
        for dem_replaced_prev in carrier_final_dem.values():
            dem_replaced = dem_replaced.drop(dem_replaced_prev.index, errors="ignore")
        carrier_final_dem[carrier] = dem_replaced

    for carrier, df in carrier_final_dem.items():
        carrier_final_dem[carrier] = (
            df["final"].to_xarray().assign_coords(carrier_name=carrier)
        )

    final_dem = xr.concat(carrier_final_dem.values(), dim="carrier_name")
    final_dem = final_dem.sum(["section", "subsection"])

    final_dem = jrc.standardize(final_dem, "twh")

    return final_dem


if __name__ == "__main__":
    get_generic_demand(
        non_generic_categories=snakemake.params.non_generic_categories,
        generic_config=snakemake.params.generic_config,
        path_energy_balances=snakemake.input.path_energy_balances,
        path_cat_names=snakemake.input.path_cat_names,
        path_carrier_names=snakemake.input.path_carrier_names,
        path_jrc_industry_energy=snakemake.input.path_jrc_industry_energy,
        path_jrc_industry_production=snakemake.input.path_jrc_industry_production,
        path_output=snakemake.output.path_output,
    )
