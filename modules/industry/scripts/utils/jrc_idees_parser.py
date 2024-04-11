import pandas as pd

MAX_YEAR = 2015


def get_auxiliary_electricity_consumption(
    process: str, jrc_energy_df: pd.DataFrame, jrc_prod_df: pd.DataFrame
) -> pd.DataFrame:
    """Get energy consumption for key industrial auxiliary processes."""
    auxiliaries = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    consumption = (
        jrc_energy_df.xs(("consumption", process))
        .loc[auxiliaries]
        .sum(level=["country_code", "unit", "cat_name"])
    )
    specific_consumption = consumption.div(jrc_prod_df.loc[process].droplevel("unit"))
    specific_consumption.index = specific_consumption.index.set_levels(
        ["ktoe/kt"], level="unit"
    )
    return specific_consumption.fillna(0)


def get_specific_electricity_consumption(
    process: str,
    subprocess: str,
    jrc_energy_df: pd.DataFrame,
    jrc_prod_df: pd.DataFrame,
) -> pd.DataFrame:
    """Get specific electricity consumption for a given JRC process -> sub-process."""
    consumption = jrc_energy_df.xs(("consumption", process, subprocess))
    demand = jrc_energy_df.xs(("demand", process, subprocess))
    specific_demand = demand.sum(level=["country_code"]).div(
        jrc_prod_df.loc[process].droplevel("unit")
    )
    efficiency = demand.div(consumption)
    electrical_efficiency = (
        efficiency.where(efficiency > 0)
        .xs("Electricity", level="carrier_name")
        .T.fillna(  # have to reorient the array thanks to a NotImplementedError
            efficiency.xs("Electricity", level="carrier_name").mean(axis=1)
        )
        .T.fillna(  # have to reorient the array thanks to a NotImplementedError
            efficiency.xs("Electricity", level="carrier_name").mean()
        )
    )

    specific_consumption = specific_demand.div(electrical_efficiency).rename(
        index={"ktoe": "ktoe/kt"}
    )
    assert (
        (
            specific_consumption.fillna(-1).droplevel("unit")
            >= specific_demand.fillna(-1)
        )
        .all()
        .all()
    )

    return specific_consumption.fillna(0)


def get_carrier_demand(
    carrier: str, all_demand_df: pd.DataFrame, jrc_energy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Get demand for a specific carrier, assuming all end use demand that could consume
    that carrier are completely met by that carrier.
    """
    energy = jrc_energy_df.xs(carrier, level="carrier_name")
    energy_efficiency = energy.xs("demand").div(energy.xs("consumption"))
    # Fill NaNs (where there is demand, but no consumption in that country)
    # with the average efficiency a. from the country, b. from all countries
    energy_efficiency = energy_efficiency.fillna(energy_efficiency.mean())

    return all_demand_df.reindex(energy_efficiency.index).div(energy_efficiency)
