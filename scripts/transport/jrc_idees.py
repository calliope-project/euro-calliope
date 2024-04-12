from pathlib import Path

import numpy as np
import pandas as pd
from eurocalliopelib import utils
from styleframe import StyleFrame

idx = pd.IndexSlice

DATASET_PARAMS = {
    "road-energy": {
        "sheet_name": "TrRoad_ene",
        "idx_start_str": "Total energy consumption",
        "idx_end_str": "Indicators",
        "unit": "ktoe",
    },
    "road-distance": {
        "sheet_name": "TrRoad_act",
        "idx_start_str": "Vehicle-km driven",
        "idx_end_str": "Stock of vehicles",
        "unit": "mio_km",
    },
    "road-vehicles": {
        "sheet_name": "TrRoad_act",
        "idx_start_str": "Stock of vehicles - in use",
        "idx_end_str": "New vehicle-registrations",
        "unit": "N. vehicles",
    },
}

ROAD_CARRIERS = {
    "Gasoline engine": "petrol",
    "Diesel oil engine": "diesel",
    "Natural gas engine": "natural_gas",
    "LPG engine": "lpg",
    "Battery electric vehicles": "electricity",
    "Plug-in hybrid electric": "petrol",
}


def process_jrc_transport_data(
    path_to_data: str,
    dataset: object,
    out_path: str,
    vehicle_type_names: dict[str, str],
) -> None:
    processed_data = pd.concat([
        read_transport_excel(path, **DATASET_PARAMS[dataset])
        for path in Path(path_to_data).glob("*.xlsx")
    ])
    if DATASET_PARAMS[dataset]["unit"] == "ktoe":
        processed_data = processed_data.apply(utils.ktoe_to_twh)
    processed_data = (
        processed_data.reset_index(level="country_code")
        .assign(country_code=lambda df: df.country_code.map(utils.convert_country_code))
        .set_index("country_code", append=True)
        .stack("year")
        .rename("value")
        .rename(index=vehicle_type_names, level="vehicle_type")
        .to_csv(out_path)
    )


def read_transport_excel(
    path: Path, sheet_name: str, idx_start_str: str, idx_end_str: str, **kwargs: object
) -> pd.DataFrame:
    xls = pd.ExcelFile(path)
    style_df = StyleFrame.read_excel(xls, read_style=True, sheet_name=sheet_name)
    df = pd.read_excel(xls, sheet_name=sheet_name)
    column_names = str(style_df.data_df.columns[0])
    # We have manually identified the section of data which is of use to us,
    # given by idx_start_str and idx_end_str.
    idx_start = int(style_df[style_df[column_names].str.find(idx_start_str) > -1][0])
    idx_end = int(style_df[style_df[column_names].str.find(idx_end_str) > -1][0])
    df = df.assign(indent=style_df[column_names].style.indent.astype(int)).loc[
        idx_start:idx_end
    ]

    total_to_check = df.iloc[0]
    # The indent of the strings in the first column of data indicates their hierarchy in a multi-level index.
    # Two levels of the hierarchy are identified here and ffill() is used to match all relevant rows
    # to the top-level name.
    df["section"] = df.where(df.indent == 1).iloc[:, 0].ffill()
    df["vehicle_type"] = df.where(df.indent == 2).iloc[:, 0].ffill()

    if sheet_name == "TrRoad_act":
        df = process_road_vehicles(df, column_names)
    elif sheet_name == "TrRoad_ene":
        df = process_road_energy(df, column_names)
    df = df.assign(country_code=column_names.split(" - ")[0]).set_index(
        "country_code", append=True
    )
    df.columns = df.columns.astype(int).rename("year")

    # After all this hardcoded cleanup, make sure numbers match up
    assert np.allclose(df.sum(), total_to_check.loc[df.columns].astype(float))

    return df


def process_road_vehicles(df: pd.DataFrame, column_names: str) -> pd.DataFrame:
    # The indent of the strings in the first column of data indicates their hierarchy in a multi-level index.
    # The vehicle subtype is identified here
    df["vehicle_subtype"] = df.where(df.indent == 3).iloc[:, 0]
    # 2-wheelers are powered by fuel oil in the dataset.
    # All useful information is either when the index column string is indented 3 times,
    # or when the vehicle type is a 2-wheeler. One of the many ways in which this dataset is a pain.
    df = df.where((df.indent == 3) | (df.vehicle_type == "Powered 2-wheelers")).dropna(
        how="all"
    )
    df.loc[df.vehicle_type == "Powered 2-wheelers", "vehicle_subtype"] = (
        "Gasoline engine"
    )
    return df.set_index(["section", "vehicle_type", "vehicle_subtype"]).drop(
        [column_names, "indent"], axis=1
    )


def process_road_energy(df: pd.DataFrame, column_names: str) -> pd.DataFrame:
    # The indent of the strings in the first column of data indicates their hierarchy in a multi-level index.
    # The vehicle subtype is identified here and ffill() is used to match all relevant rows to this subtype.
    df["vehicle_subtype"] = df.where(df.indent == 3).iloc[:, 0].ffill()
    # Remove bracketed information from the vehicle type and subtype
    df["vehicle_type"] = df["vehicle_type"].str.split("(", expand=True)[0].str.strip()
    df["vehicle_subtype"] = (
        df["vehicle_subtype"].str.split("(", expand=True)[0].str.strip()
    )
    # Powered 2-wheelers are gasoline engine only (this is implicit when looking at the Excel sheet directly)
    df.loc[df.vehicle_type == "Powered 2-wheelers", "vehicle_subtype"] = (
        "Gasoline engine"
    )
    df["carrier"] = df.where(df.indent == 4).iloc[:, 0]
    df["carrier"] = df["carrier"].str.replace("of which ", "")
    # Powered 2-wheelers use petrol, some of which is biofuels (we deal with the 'of which' part later)
    df.loc[(df.vehicle_type == "Powered 2-wheelers") & (df.indent == 2), "carrier"] = (
        "petrol"
    )
    df.loc[(df.vehicle_type == "Powered 2-wheelers") & (df.indent == 3), "carrier"] = (
        "biofuels"
    )
    # both domestic and international freight uses diesel in the dataset
    # (this is implicit when looking at the Excel sheet directly)
    df.loc[(df.vehicle_subtype == "Domestic") & (df.indent == 3), "carrier"] = "diesel"
    df.loc[(df.vehicle_subtype == "International") & (df.indent == 3), "carrier"] = (
        "diesel"
    )
    # All other vehicle types mention the drive-train directly, so we translate that to energy carrier here
    df["carrier"] = df["carrier"].fillna(df.vehicle_subtype.replace(ROAD_CARRIERS))

    df = (
        df.where((df.indent > 2) | (df.vehicle_type == "Powered 2-wheelers"))
        .dropna()
        .set_index(["section", "vehicle_type", "vehicle_subtype", "carrier"])
        .drop([column_names, "indent"], axis=1)
    )

    df = remove_of_which(df, "diesel", "biofuels")
    df = remove_of_which(df, "petrol", "biofuels")
    df = remove_of_which(df, "petrol", "electricity")
    df = remove_of_which(df, "natural_gas", "biogas")

    return df


def remove_of_which(
    df: pd.DataFrame, main_carrier: str, of_which_carrier: str
) -> pd.DataFrame:
    """
    Subfuels (e.g. biodiesel) are given as 'of which ...', meaning the main fuel consumption
    includes those fuels (diesel is actually diesel + biofuels). We rectify that here
    """
    updated_carrier = (
        (
            df.xs(main_carrier, level="carrier")
            - df.xs(of_which_carrier, level="carrier")
        )
        .dropna()
        .assign(carrier=main_carrier)
        .set_index("carrier", append=True)
    )
    df.loc[updated_carrier.index] = updated_carrier
    return df


if __name__ == "__main__":
    process_jrc_transport_data(
        path_to_data=snakemake.input.data,
        dataset=snakemake.wildcards.dataset,
        out_path=snakemake.output[0],
        vehicle_type_names={
            v: k for k, v in snakemake.params.vehicle_type_names.items()
        },
    )
