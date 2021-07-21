from pathlib import Path

import numpy as np
import pandas as pd
from styleframe import StyleFrame

from eurocalliopelib import utils

idx = pd.IndexSlice

DATASET_PARAMS = {
    "road-energy": {
        "sheet_name": "TrRoad_ene",
        "idx_start_str": "Total energy consumption",
        "idx_end_str": "Indicators",
        "unit": "ktoe"
    },
    "road-distance": {
        "sheet_name": "TrRoad_act",
        "idx_start_str": "Vehicle-km driven",
        "idx_end_str": "Stock of vehicles",
        "unit": "mio_km"
    },
    "road-vehicles": {
        "sheet_name": "TrRoad_act",
        "idx_start_str": "Stock of vehicles - in use",
        "idx_end_str": "New vehicle-registrations",
        "unit": "N. vehicles"
    },
    "rail-energy": {
        "sheet_name": "TrRail_ene",
        "idx_start_str": "Total energy consumption",
        "idx_end_str": "Indicators",
        "unit": "ktoe"
    },
    "rail-distance": {
        "sheet_name": "TrRail_act",
        "idx_start_str": "Vehicle-km (mio km)",
        "idx_end_str": "Stock of vehicles",
        "unit": "mio_km"
    },
}

ROAD_CARRIERS = {
    'Gasoline engine': 'petrol',
    'Diesel oil engine': 'diesel',
    'Natural gas engine': 'natural_gas',
    'LPG engine': 'lpg',
    'Battery electric vehicles': 'electricity',
    'Plug-in hybrid electric': 'petrol'
}

RAIL_CARRIERS = {
    'Diesel oil (incl. biofuels)': 'diesel',
    'Electric': 'electricity',
    'Diesel oil': 'diesel'
}


def process_jrc_transport_data(data_dir, dataset, out_path):
    data_filepaths = list(Path(data_dir).glob("*.xlsx"))
    processed_data = pd.concat([
        read_transport_excel(file, **DATASET_PARAMS[dataset])
        for file in data_filepaths
    ])
    if DATASET_PARAMS[dataset]["unit"] == "ktoe":
        processed_data = processed_data.apply(utils.ktoe_to_twh)
        processed_data.index = processed_data.index.set_levels(['twh'], level='unit')

    processed_data.stack('year').to_csv(out_path)


def read_transport_excel(file, sheet_name, idx_start_str, idx_end_str, unit):
    style_df = StyleFrame.read_excel(file, read_style=True, sheet_name=sheet_name)
    df = pd.read_excel(file, sheet_name=sheet_name)
    column_names = str(style_df.data_df.columns[0])
    idx_start = int(style_df[style_df[column_names].str.find(idx_start_str) > -1][0])
    idx_end = int(style_df[style_df[column_names].str.find(idx_end_str) > -1][0])
    df = df.assign(indent=style_df[column_names].style.indent.astype(int)).loc[idx_start:idx_end]

    total_to_check = df.iloc[0]
    df['section'] = df.where(df.indent == 1).iloc[:, 0].ffill()
    df['vehicle_type'] = df.where(df.indent == 2).iloc[:, 0].ffill()

    if sheet_name == 'TrRoad_act':
        df = process_road_vehicles(df, column_names)
    elif sheet_name == 'TrRoad_ene':
        df = process_road_energy(df, column_names)
    elif sheet_name == 'TrRail_ene' or sheet_name == 'TrRail_act':
        df = process_rail(df, column_names)

    df = (
        df
        .assign(country_code=column_names.split(' - ')[0], unit=unit)
        .set_index(['country_code', 'unit'], append=True)
    )
    df.columns = df.columns.astype(int).rename('year')

    # After all this hardcoded cleanup, make sure numbers match up
    assert np.allclose(
        df.sum(),
        total_to_check.loc[df.columns].astype(float)
    )

    return df


def process_rail(df, column_names):
    df['carrier'] = df.where(df.indent == 3).iloc[:, 0]
    # ASSUME: All metro/tram/high speed rail is electrically powered
    df.loc[df.vehicle_type == 'Metro and tram, urban light rail', 'carrier'] = 'electricity'
    df.loc[df.vehicle_type == 'High speed passenger trains', 'carrier'] = 'electricity'

    df['carrier'] = df['carrier'].replace(RAIL_CARRIERS).fillna(df.vehicle_type.replace(RAIL_CARRIERS))
    df.loc[df.section == 'Freight transport', 'vehicle_type'] = 'Freight'

    return (
        df
        .where((df.indent > 1) & (df.carrier.str.find('Conventional') == -1))
        .dropna()
        .set_index(['section', 'vehicle_type', 'carrier'])
        .drop([column_names, 'indent'], axis=1)
    )


def process_road_vehicles(df, column_names):
    df['vehicle_subtype'] = df.where(df.indent == 3).iloc[:, 0]
    # ASSUME: 2-wheelers are powered by fuel oil
    df = df.where(
        (df.indent == 3) | (df.vehicle_type == 'Powered 2-wheelers')
    ).dropna(how='all')
    df.loc[df.vehicle_type == 'Powered 2-wheelers', 'vehicle_subtype'] = 'Gasoline engine'
    return (
        df
        .set_index(['section', 'vehicle_type', 'vehicle_subtype'])
        .drop([column_names, 'indent'], axis=1)
    )


def process_road_energy(df, column_names):
    df['vehicle_subtype'] = df.where(df.indent == 3).iloc[:, 0].ffill()
    df['vehicle_type'] = df['vehicle_type'].str.split('(', expand=True)[0].str.strip()
    df['vehicle_subtype'] = df['vehicle_subtype'].str.split('(', expand=True)[0].str.strip()
    df.loc[df.vehicle_type == 'Powered 2-wheelers', 'vehicle_subtype'] = 'Gasoline engine'
    df['carrier'] = df.where(df.indent == 4).iloc[:, 0]
    df['carrier'] = df['carrier'].str.replace('of which ', '')
    df.loc[(df.vehicle_type == 'Powered 2-wheelers') & (df.indent == 2), 'carrier'] = 'petrol'
    df.loc[(df.vehicle_type == 'Powered 2-wheelers') & (df.indent == 3), 'carrier'] = 'biofuels'
    df.loc[(df.vehicle_subtype == 'Domestic') & (df.indent == 3), 'carrier'] = 'diesel'
    df.loc[(df.vehicle_subtype == 'International') & (df.indent == 3), 'carrier'] = 'diesel'
    df['carrier'] = df['carrier'].fillna(df.vehicle_subtype.replace(ROAD_CARRIERS))
    df = (
        df
        .where((df.indent > 2) | (df.vehicle_type == 'Powered 2-wheelers'))
        .dropna()
        .set_index(['section', 'vehicle_type', 'vehicle_subtype', 'carrier'])
        .drop([column_names, 'indent'], axis=1)
    )

    df = remove_of_which(df, 'diesel', 'biofuels')
    df = remove_of_which(df, 'petrol', 'biofuels')
    df = remove_of_which(df, 'petrol', 'electricity')
    df = remove_of_which(df, 'natural_gas', 'biogas')

    return df


def remove_of_which(df, main_carrier, of_which_carrier):
    """
    Subfuels (e.g. biodiesel) are given as 'of which ...', meaning the main fuel consumption
    includes those fuels (diesel is actually diesel + biofuels). We rectify that here
    """
    updated_carrier = (
        (df.xs(main_carrier, level='carrier') - df.xs(of_which_carrier, level='carrier'))
        .dropna()
        .assign(carrier=main_carrier)
        .set_index('carrier', append=True)
    )
    df.loc[updated_carrier.index] = updated_carrier
    return df


if __name__ == "__main__":
    process_jrc_transport_data(
        data_dir=snakemake.input.unprocessed_data,
        dataset=snakemake.wildcards.dataset,
        out_path=snakemake.output[0]
    )
