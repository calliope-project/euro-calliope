from typing import Optional

import numpy as np
import pandas as pd
from eurocalliopelib import utils

CH_ENERGY_CARRIER_TRANSLATION = {
    "Heizöl": "oil",
    "Erdgas": "gas",
    "El. Widerstandsheizungen": "direct_electric",
    "El. Wärmepumpen 1)": "heat_pump",
    "El. Ohm'sche Anlagen": "direct_electric",
    "El. Wärmepumpen": "heat_pump",
    "Elektrizität": "electricity",
    "Holz": "biofuel",
    "Kohle": "solid_fossil",
    "Fernwärme": "heat",
    "Umweltwärme": "ambient_heat",
    "Solar": "solar_thermal",
}

CH_HH_END_USE_TRANSLATION = {
    "Raumwärme": "space_heat",
    "Warmwasser": "water_heat",
    "Prozesswärme": "cooking",
    "Beleuchtung": "end_use_electricity",
    "Klima, Lüftung, HT": "end_use_electricity",
    "I&K, Unterhaltung": "end_use_electricity",
    "Antriebe, Prozesse": "end_use_electricity",
    "sonstige": "end_use_electricity",
}


def ch_building_energy_consumption_dataset(
    path_to_ch_end_use_excel: str, path_to_output: str
):
    household_consumption = ch_hh_consumption(path_to_ch_end_use_excel)
    commercial_consumption = ch_non_hh_consumption(
        path_to_ch_end_use_excel, household_consumption
    )

    building_energy_consumption = utils.merge_da(
        [
            household_consumption.expand_dims(cat_name=["household"]),
            commercial_consumption.expand_dims(cat_name=["commercial"]),
        ],
        merged_da_name="building_heat_demand",
    )

    final_da = (
        utils.pj_to_twh(building_energy_consumption)
        .expand_dims(country_code=["CHE"])
        .assign_attrs(unit="twh")
    )
    final_da.to_netcdf(path_to_output)


def ch_non_hh_consumption(path_to_ch_end_use_excel: str, household_consumption: str):
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """

    commercial_sheet_kwargs = {
        "translation": CH_HH_END_USE_TRANSLATION,
        "index_name": "end_use",
    }
    commercial_fuel_consumption = get_ch_sheet(
        path_to_ch_end_use_excel, "Tabelle 25", skipfooter=4, **commercial_sheet_kwargs
    )
    # Quirk of the excel is that there is no space in this sheet name
    commercial_electricity_consumption = get_ch_sheet(
        path_to_ch_end_use_excel, "Tabelle26", skipfooter=4, **commercial_sheet_kwargs
    )
    commercial_electricity_consumption = commercial_electricity_consumption.expand_dims(
        carrier_name=["electricity"]
    )

    # ASSUME: commercial fuel consumption is assigned to end uses (space heat, cooking, water heating)
    # with the same ratio as in non-commercial households (since we do not have this data).

    relevant_household_consumption = household_consumption.sel(
        end_use=commercial_fuel_consumption.end_use
    ).dropna("carrier_name")

    household_consumption_shares = (
        relevant_household_consumption
        / relevant_household_consumption.sum("carrier_name")
    )
    commercial_fuel_consumption_assigned_to_end_uses = (
        household_consumption_shares * commercial_fuel_consumption
    )

    assert np.allclose(  # ensure we haven't lost some energy demand in the process
        commercial_fuel_consumption_assigned_to_end_uses.sum("carrier_name"),
        commercial_fuel_consumption,
    )
    commercial_energy_consumption = utils.merge_da([
        commercial_electricity_consumption,
        commercial_fuel_consumption_assigned_to_end_uses,
    ])

    return commercial_energy_consumption


def ch_hh_consumption(path_to_ch_end_use_excel: str):
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """
    household_sheet_kwargs = {
        "translation": CH_ENERGY_CARRIER_TRANSLATION,
        "index_name": "carrier_name",
    }
    household_space_heat_consumption = get_ch_sheet(
        path_to_ch_end_use_excel, "Tabelle 18", skipfooter=8, **household_sheet_kwargs
    )
    household_hot_water_consumption = get_ch_sheet(
        path_to_ch_end_use_excel, "Tabelle 20", skipfooter=5, **household_sheet_kwargs
    )
    # Quirk of the excel is that there is no space in this sheet name
    household_cooking_consumption = get_ch_sheet(
        path_to_ch_end_use_excel, "Tabelle21", skipfooter=4, **household_sheet_kwargs
    )

    household_energy_consumption = utils.merge_da([
        household_space_heat_consumption.expand_dims(end_use=["space_heat"]),
        household_hot_water_consumption.expand_dims(end_use=["water_heat"]),
        household_cooking_consumption.expand_dims(end_use=["cooking"]),
    ])

    return household_energy_consumption


def get_ch_sheet(
    path_to_excel: str,
    sheet: str,
    skipfooter: int,
    index_name: str,
    translation: Optional[dict] = None,
):
    df = pd.read_excel(
        path_to_excel, sheet_name=sheet, skiprows=9, skipfooter=skipfooter, index_col=1
    )
    to_drop = [col for col in df.columns if col.startswith(("Unnamed", "Δ"))]
    df = df.drop(to_drop, axis=1)
    df.index = df.index.str.strip()
    df.columns = df.columns.astype(int)

    if translation is not None:
        df = df.groupby(translation).sum()

    df = df.rename_axis(index=index_name, columns="year")

    return df.stack().to_xarray()


if __name__ == "__main__":
    ch_building_energy_consumption_dataset(
        path_to_ch_end_use_excel=snakemake.input.ch_end_use_excel,
        path_to_output=snakemake.output[0],
    )
