from typing import Optional

import pandas as pd
import xarray as xr
from utils import formatting
from utils import jrc_idees_parser as jrc

# TODO: this should be defined externally via a .csv file.
JRC_TO_CALLIOPE = {
    "Electricity": "electricity",
    "Natural gas (incl. biogas)": "methane",
    "Diesel oil (incl. biofuels)": "diesel",
    "Low enthalpy heat": "space_heat",
}


def get_other_demand(
    config_params: dict,
    path_energy_balances: str,
    path_cat_names: str,
    path_carrier_names: str,
    path_jrc_industry_energy: str,
    path_jrc_industry_production: str,
    path_output: Optional[str] = None,
) -> pd.DataFrame:
    """Execute the default data processing pipeline all non-specific industries.

    Args:
        config_params (dict): all industry configuration parameters.
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
        path_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True
    )
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    jrc_energy = xr.open_dataset(path_jrc_industry_energy)
    jrc_prod = xr.open_dataset(path_jrc_industry_production)

    # TODO: fix naming convention forced by the JRC module.
    jrc_energy = jrc_energy.rename({"jrc-idees-industry-twh": "value"})
    jrc_prod = jrc_prod.rename({"jrc-idees-industry-twh": "value"})

    # Remove data from all specifically processed industries
    specific_industries = config_params["specific-industries"]
    cat_names_df = cat_names_df[~cat_names_df["jrc_idees"].isin(specific_industries)]
    jrc_energy = jrc_energy.drop_sel(cat_name=specific_industries)
    jrc_prod = jrc_prod.drop_sel(cat_name=specific_industries)

    # Process data:
    # Extract useful dem. -> remove useful dem. from rest -> extract final dem.
    selected_useful = config_params["other"]["useful-demands"]
    other_useful_demand = extract_aggregated_useful_demand(jrc_energy, selected_useful)
    other_useful_demand = other_useful_demand.rename({"subsection": "carrier_name"})

    selected_final = config_params["other"]["final-energy-carriers"]
    final_method = config_params["other"]["final-energy-method"]
    jrc_energy_clean = jrc_energy.drop_sel(subsection=selected_useful)
    if final_method == "priority":
        other_final_demand = extract_final_demand_by_priority(
            jrc_energy_clean, selected_final
        )
    else:
        raise ValueError(f"Unsupported final energy method: {final_method}.")

    # Combine and fill missing countries
    other_demand = xr.concat(
        [other_useful_demand, other_final_demand], dim="carrier_name"
    )
    other_demand = formatting.fill_missing_countries_years(
        energy_balances_df, cat_names_df, carrier_names_df, other_demand
    )

    # Fix the naming
    for carrier in JRC_TO_CALLIOPE:
        other_demand["carrier_name"] = other_demand["carrier_name"].where(
            lambda x, i=carrier: x != i, JRC_TO_CALLIOPE[carrier]
        )
    other_demand = jrc.ensure_standard_coordinates(other_demand)
    other_demand["value"].attrs["units"] = "twh"

    if path_output:
        other_demand.to_netcdf(path_output)

    return other_demand


def extract_aggregated_useful_demand(
    jrc_energy: xr.Dataset, subsections: list[str]
) -> xr.Dataset:
    """Gather useful demand of specific subsections and aggregate them.

    Args:
        jrc_energy (xr.Dataset): JRC energy dataset.
        subsections (list[str]): subsections to take.

    Returns:
        xr.Dataset: dataset with [subsection, cat_name, year, country_code] dimensions.
    """
    useful_dem_total = jrc_energy.sel(energy="demand").sum("carrier_name")

    useful_dem_extracted = useful_dem_total.sel(subsection=subsections)
    useful_dem_extracted = useful_dem_extracted.sum("section").drop("energy")

    return useful_dem_extracted


def extract_final_demand_by_priority(
    jrc_energy: xr.Dataset, carrier_priority: list[str]
):
    """Gather final demand of all sectors giving priority to certain carriers.

    For the given carriers: drop overlapping consumption so that any demand
    is met by carriers with the given priority.
    E.g., if carrier priority is [Electricity, Natural gas, Diesel] then:
    - Electricity: if met exclusively or otherwise, it's final electrical demand.
    - Natural gas: if met exclusively or otherwise, EXCEPT for overlapping cases with the above.
    - Diesel: if met exclusively or otherwise, EXCEPT for overlapping cases with all the above.

    Args:
        jrc_energy (xr.Dataset): JRC energy dataset.
        carrier_priority (list[str]): carriers to take in order of priority.

    Returns:
        xr.Dataset: dataset filled with demands for the given carriers.
    """
    useful_dem = jrc_energy.sel(energy="demand").sum("carrier_name")

    final_dem_carrier = {}
    for carrier in carrier_priority:
        dem_current = jrc.get_carrier_final_demand(carrier, useful_dem, jrc_energy)
        dem_current = dem_current.to_dataframe().dropna()
        for dem_prev in final_dem_carrier.values():
            dem_current = dem_current.drop(dem_prev.index, errors="ignore")
        final_dem_carrier[carrier] = dem_current

    for carrier, df in final_dem_carrier.items():
        final_dem_carrier[carrier] = df.to_xarray().assign_coords(carrier_name=carrier)

    final_dem = xr.concat(final_dem_carrier.values(), dim="carrier_name")
    final_dem = final_dem.sum(["section", "subsection"])

    return final_dem


if __name__ == "__main__":
    get_other_demand(
        config_params=snakemake.params.config_params,
        path_energy_balances=snakemake.input.path_energy_balances,
        path_cat_names=snakemake.input.path_cat_names,
        path_carrier_names=snakemake.input.path_carrier_names,
        path_jrc_industry_energy=snakemake.input.path_jrc_industry_energy,
        path_jrc_industry_production=snakemake.input.path_jrc_industry_production,
        path_output=snakemake.output.path_output,
    )
