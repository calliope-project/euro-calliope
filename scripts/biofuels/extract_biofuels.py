import pandas as pd

from eurocalliopelib import utils

SHEET_NAME_POTENTIALS = "ENER - NUTS0 EnergyCom"
SHEET_NAME_COSTS = "COST - NUTS0 EnergyCom"

COL_NAME_SCENARIO = "Scenario"
COL_NAME_FEEDSTOCK = "Energy Commodity"
COL_NAME_COUNTRY = "NUTS0"
COL_NAME_YEAR = "Year"
COL_NAME_VALUE = "Value"
COL_NAME_VALUE2 = "NUTS0 Energy Commodity Cost "

SCENARIO_NAME_MAP = {
    "ENS_Low": "low",
    "ENS_Med": "medium",
    "ENS_High": "high"
}

FEEDSTOCK_NAME_MAP = {
    "MINBIOFRSR1": "forestry-energy-residues",
    "MINBIOFRSR1a": "landscape-care-residues",
    "MINBIOGAS1": "manure",
    "MINBIOMUN1": "municipal-waste",
    "MINBIOAGRW1": "primary-agricultural-residues",
    "MINBIOWOOa": "roundwood-chips",
    "MINBIOWOO": "roundwood-fuelwood",
    "MINBIOWOOW1a": "secondary-forestry-residues-sawdust",
    "MINBIOWOOW1": "secondary-forestry-residues-woodchips",
    "MINBIOSLU1": "sludge"
}

COL_NAME_MAP = {
    COL_NAME_COUNTRY: "country_code",
    COL_NAME_FEEDSTOCK: "feedstock",
    COL_NAME_SCENARIO: "scenario",
    COL_NAME_YEAR: "year",
    COL_NAME_VALUE: "value",
    COL_NAME_VALUE2: "value"
}


def extract_biofuels(path_to_raw_data, paths_to_outputs):
    potentials = read_filter_and_map_names(path_to_raw_data, SHEET_NAME_POTENTIALS)
    potentials.to_csv(paths_to_outputs.potentials, index=True, header=True)
    costs = read_filter_and_map_names(path_to_raw_data, SHEET_NAME_COSTS)
    costs.to_csv(paths_to_outputs.costs, index=True, header=True)


def read_filter_and_map_names(path_to_raw_data, sheet_name):
    # ASSUME no energy crop feedstock available.
    # Ignores scenarios other than high/med/low.
    # Removes "KS" NUTS0 from costs, as this is an unknown abbreviation.
    # Renames scenarios and feedstocks.
    # Renames countries to ISO3.
    # Renames column names
    # Sets 1 missing value (BIH med potential 2010) to 0.
    # Sets 1 missing value (NOR, high cost 2030) to 0.
    df = pd.read_excel(path_to_raw_data, sheet_name=sheet_name, na_values="")
    df = df[df[COL_NAME_FEEDSTOCK].isin(FEEDSTOCK_NAME_MAP.keys())]
    df = df[df[COL_NAME_SCENARIO].isin(SCENARIO_NAME_MAP.keys())]
    df = df.where(df[COL_NAME_COUNTRY] != "KS").dropna(subset=[COL_NAME_COUNTRY])
    df[COL_NAME_SCENARIO] = df[COL_NAME_SCENARIO].map(SCENARIO_NAME_MAP)
    df[COL_NAME_FEEDSTOCK] = df[COL_NAME_FEEDSTOCK].map(FEEDSTOCK_NAME_MAP)
    df[COL_NAME_COUNTRY] = df[COL_NAME_COUNTRY].map(utils.eu_country_code_to_iso3)
    df.rename(columns=COL_NAME_MAP, inplace=True)
    df.value = pd.to_numeric(df.value, errors='coerce')
    assert df.value.isna().sum() <= 1
    df.value.fillna(0, inplace=True)
    assert df.value.isna().sum() == 0
    return df


if __name__ == "__main__":
    extract_biofuels(
        path_to_raw_data=snakemake.input.potentials_and_costs,
        paths_to_outputs=snakemake.output
    )
