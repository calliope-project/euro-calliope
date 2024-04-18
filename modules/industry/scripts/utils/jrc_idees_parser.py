import pandas as pd
import xarray as xr

MAX_YEAR = 2015


def get_specific_final_demand_electric_auxiliary(
    section: str, jrc_energy: xr.Dataset, jrc_prod: xr.Dataset
) -> pd.DataFrame:
    """Get auxiliary electricity consumption for auxiliary processes."""
    # Extract relevant sector/subsector data
    auxiliaries = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    final_demand = jrc_energy.sel(
        energy="consumption",
        carrier_name="Electricity",
        section=section,
        subsection=auxiliaries,
    )
    production = jrc_prod.sel(produced_material=section)

    # Calculate specific final demand (twh/kt) of auxiliary processes
    final_demand = final_demand.sum(dim="subsection")
    specific_final_demand = final_demand / production

    specific_final_demand = specific_final_demand["value"].assign_attrs({
        "units": "twh/kt"
    })

    return specific_final_demand.fillna(0)


def get_specific_final_demand_electric(
    section: str,
    subsection: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.Dataset,
) -> pd.DataFrame:
    """Get specific electricity consumption for a given JRC section -> subsection."""
    # Extract relevant section and subsection data.
    final_demand = jrc_energy.sel(
        energy="consumption", section=section, subsection=subsection
    )
    useful_demand = jrc_energy.sel(
        energy="demand", section=section, subsection=subsection
    )
    production = jrc_prod.sel(produced_material=section)

    # Calculate electrical efficiency and fill empty values
    # First by country avg. over all years, then by year avg. over all countries.
    eff = useful_demand / final_demand
    eff_electric = eff.where(eff > 0).sel(carrier_name="Electricity")
    eff_electric = eff_electric.fillna(eff_electric.mean(dim="year"))
    eff_electric = eff_electric.fillna(eff_electric.mean(dim="country_code"))

    # Get the useful energy demand per production method (e.g., twh/kt_steel)
    specific_useful_demand = useful_demand.sum(dim=["carrier_name"]) / production
    # Use the efficiency to determine final electrical demand.
    specific_final_demand = specific_useful_demand / eff_electric
    specific_final_demand = specific_final_demand["value"].assign_attrs({
        "units": "twh/kt"
    })

    assert specific_final_demand >= specific_useful_demand, "Creating energy!"

    return specific_final_demand.fillna(0)


# TODO: fix me!
def get_carrier_demand(
    carrier: str, all_demand_df: pd.DataFrame, jrc_energy: xr.Dataset
) -> pd.DataFrame:
    """
    Get demand for a specific carrier, assuming all end use demand that could consume
    that carrier are completely met by that carrier.
    """
    energy = jrc_energy.xs(carrier, level="carrier_name")
    energy_efficiency = energy.xs("demand").div(energy.xs("consumption"))
    # Fill NaNs (where there is demand, but no consumption in that country)
    # with the average efficiency a. from the country, b. from all countries
    energy_efficiency = energy_efficiency.fillna(energy_efficiency.mean())

    return all_demand_df.reindex(energy_efficiency.index).div(energy_efficiency)
