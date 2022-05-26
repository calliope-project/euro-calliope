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

COL_NAME_MAP = {
    COL_NAME_COUNTRY: "country_code",
    COL_NAME_FEEDSTOCK: "feedstock",
    COL_NAME_SCENARIO: "scenario",
    COL_NAME_YEAR: "year",
    COL_NAME_VALUE: "value",
    COL_NAME_VALUE2: "value"
}

MISSING_COST_VALUE = (2030, "high", "NOR", "municipal-waste")
REPLACEMENT_COST_VALUE = (2030, "medium", "NOR", "municipal-waste")

MISSING_POTENTIAL_VALUE = (2010, "medium", "BIH", "sludge")
REPLACEMENT_POTENTIAL_VALUE = (2010, "high", "BIH", "sludge") # 0 PJ


def extract_biofuels(path_to_raw_data, feedstocks, paths_to_outputs):
    potentials = read_and_filter_and_map_names(path_to_raw_data, SHEET_NAME_POTENTIALS, feedstocks)
    potentials = replace(potentials, MISSING_POTENTIAL_VALUE, REPLACEMENT_POTENTIAL_VALUE)
    assert potentials["value"].isna().sum() == 0
    potentials.to_csv(paths_to_outputs.potentials, index=True, header=True)

    costs = read_and_filter_and_map_names(path_to_raw_data, SHEET_NAME_COSTS, feedstocks)
    costs = replace(costs, MISSING_COST_VALUE, REPLACEMENT_COST_VALUE)
    assert costs["value"].isna().sum() == 0
    costs.to_csv(paths_to_outputs.costs, index=True, header=True)


def read_and_filter_and_map_names(path_to_raw_data, sheet_name, feedstocks):
    # Ignores scenarios other than high/med/low.
    # Removes "KS" NUTS0 from costs, as this is an unknown abbreviation.
    # Renames scenarios and feedstocks.
    # Renames countries to ISO3.
    # Renames column names.
    df = pd.read_excel(path_to_raw_data, sheet_name=sheet_name, na_values="")
    df = df[df[COL_NAME_FEEDSTOCK].isin(feedstocks.keys())]
    df = df[df[COL_NAME_SCENARIO].isin(SCENARIO_NAME_MAP.keys())]
    df = df.where(df[COL_NAME_COUNTRY] != "KS").dropna(subset=[COL_NAME_COUNTRY])
    df[COL_NAME_SCENARIO] = df[COL_NAME_SCENARIO].map(SCENARIO_NAME_MAP)
    df[COL_NAME_FEEDSTOCK] = df[COL_NAME_FEEDSTOCK].map(feedstocks)
    df[COL_NAME_COUNTRY] = df[COL_NAME_COUNTRY].map(utils.convert_country_code)
    df.rename(columns=COL_NAME_MAP, inplace=True)
    df["value"] = pd.to_numeric(df.value, errors='coerce')
    return df


def replace(df, index_to_replace, index_replace_from):
    # Replace one missing entry with another entry.
    df = df.set_index(["year", "scenario", "country_code", "feedstock"])
    assert pd.isna(df.loc[index_to_replace, "value"])
    df.loc[index_to_replace, "value"] = df.loc[index_replace_from, "value"]
    assert not pd.isna(df.loc[index_to_replace, "value"])
    return df.reset_index()


if __name__ == "__main__":
    extract_biofuels(
        path_to_raw_data=snakemake.input.potentials_and_costs,
        feedstocks=snakemake.params.feedstocks,
        paths_to_outputs=snakemake.output
    )
