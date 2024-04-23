from typing import Optional

import eurocalliopelib.utils as ec_utils
import pandas as pd
from utils import formatting as fmt
from utils import jrc_idees_parser as jrc


def get_other_demand(
    year_range: list,
    specific_industries: list[str],
    path_energy_balances: str,
    path_cat_names: str,
    path_carrier_names: str,
    path_jrc_energy: str,
    path_jrc_production: str,
    path_output: Optional[str] = None,
) -> pd.DataFrame:
    """Execute the default data processing pipeline all non-specific industries.

    Args:
        year_range (list): years to include and interpolate.
        specific_industries (list[str]): individually processed industries (omitted).
        path_energy_balances (str): country energy balances (usually from eurostat).
        path_cat_names (str): eurostat category mapping file.
        path_carrier_names (str): eurostat carrier name mapping file.
        path_jrc_energy (str): jrc country-specific industrial energy demand file.
        path_jrc_production (str): jrc country-specific industrial production file.
        path_output (str): location of steel demand output file.

    Returns:
        pd.DataFrame: dataframe with industrial demand per country.
    """
    # -------------------------------------------------------------------------
    # Prepare data files
    # -------------------------------------------------------------------------
    energy_balances_df = pd.read_csv(
        path_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True
    )
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    jrc_energy_df = pd.read_csv(path_jrc_energy, index_col=[0, 1, 2, 3, 4, 5, 6])
    jrc_prod_df = pd.read_csv(path_jrc_production, index_col=[0, 1, 2, 3])
    # Remove data from all specifically processed industries
    cat_names_df = cat_names_df[~cat_names_df["jrc_idees"].isin(specific_industries)]
    jrc_energy_df = jrc_energy_df.drop(specific_industries, level="cat_name")
    jrc_prod_df = jrc_prod_df.drop(specific_industries, level="cat_name")

    # -------------------------------------------------------------------------
    # Process data
    # -------------------------------------------------------------------------
    levels_to_sum = ["section", "subsection", "country_code", "cat_name", "unit"]
    demand = jrc_energy_df.xs("demand").sum(level=levels_to_sum)

    # TODO: this should be specified in the configuration file

    # if it can be met by electricity (exclusively or otherwise),
    # then it's an end-use electricity demand
    electrical_consumption = (
        jrc.get_carrier_demand("Electricity", demand, jrc_energy_df)
        .assign(carrier="electricity")
        .set_index("carrier", append=True)
    )
    # If it can only be met by natural gas (steam heating) then it's natural gas
    nat_gas_consumption = (
        jrc.get_carrier_demand("Natural gas (incl. biogas)", demand, jrc_energy_df)
        .drop(electrical_consumption.droplevel("carrier").index, errors="ignore")
        .assign(carrier="methane")
        .set_index("carrier", append=True)
    )
    # If it can only be met by diesel (backup generators) then it's diesel
    diesel_consumption = (
        jrc.get_carrier_demand("Diesel oil (incl. biofuels)", demand, jrc_energy_df)
        .drop(nat_gas_consumption.droplevel("carrier").index, errors="ignore")
        .drop(electrical_consumption.droplevel("carrier").index, errors="ignore")
        .assign(carrier="diesel")
        .set_index("carrier", append=True)
    )

    # Combine carrier files, removing space heating
    other_consumption_carrier_based = pd.concat([
        i.drop("Low enthalpy heat", level="subsection", errors="ignore").sum(
            level=["cat_name", "country_code", "unit", "carrier"]
        )
        for i in [electrical_consumption, nat_gas_consumption, diesel_consumption]
    ])
    other_consumption_use_based = (
        demand.xs("Low enthalpy heat", level="subsection")
        .sum(level=["cat_name", "country_code", "unit"])
        .assign(carrier="space_heat")
        .set_index("carrier", append=True)
    )

    all_other_consumption = pd.concat([
        other_consumption_carrier_based,
        other_consumption_use_based,
    ])

    # -------------------------------------------------------------------------
    # Format the final output
    # -------------------------------------------------------------------------
    all_other_consumption.columns = all_other_consumption.columns.astype(int).rename(
        "year"
    )
    all_other_consumption_filled = fmt.fill_missing_data(
        energy_balances_df,
        cat_names_df,
        carrier_names_df,
        all_other_consumption,
        year_range,
    )

    units = all_other_consumption_filled.index.get_level_values("unit")
    all_other_consumption_filled.loc[units == "ktoe"] = (
        all_other_consumption_filled.loc[units == "ktoe"].apply(ec_utils.ktoe_to_twh)
    )
    all_other_consumption_filled = all_other_consumption_filled.rename(
        {"ktoe": "twh"}, level="unit"
    )
    all_other_consumption_filled.index = all_other_consumption_filled.index.set_names(
        "subsector", level="cat_name"
    )
    all_other_consumption_filled = all_other_consumption_filled.stack()
    breakpoint()
    if path_output:
        all_other_consumption_filled.reorder_levels(fmt.LEVEL_ORDER).to_csv(path_output)

    return all_other_consumption_filled


if __name__ == "__main__":
    get_other_demand(
        year_range=snakemake.params.year_range,
        specific_industries=snakemake.params.specific_industries,
        path_energy_balances=snakemake.input.path_energy_balances,
        path_cat_names=snakemake.input.path_cat_names,
        path_carrier_names=snakemake.input.path_carrier_names,
        path_jrc_energy=snakemake.input.path_jrc_energy,
        path_jrc_production=snakemake.input.path_jrc_production,
        path_output=snakemake.output.path_output,
    )
