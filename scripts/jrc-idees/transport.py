import glob
import os

import numpy as np
import pandas as pd
from styleframe import StyleFrame

from eurocalliopelib import utils

idx = pd.IndexSlice

DATSET_PARAMS = {
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


def process_jrc_transport_data(data_dir, dataset, out_path):
    data_filepaths = glob.glob(os.path.join(data_dir, "*.xlsx"))
    processed_data = read_transport_excel(
        data_filepaths, **DATSET_PARAMS[dataset]
    )
    if DATSET_PARAMS[dataset]["unit"] == "ktoe":
        processed_data = processed_data.apply(utils.ktoe_to_twh)
        processed_data.index = processed_data.index.set_levels(['twh'], level='unit')

    processed_data.stack('year').to_csv(out_path)


def read_transport_excel(files, sheet_name, idx_start_str, idx_end_str, unit):
    dfs = []
    for file in files:
        style_df = StyleFrame.read_excel(file, read_style=True, sheet_name=sheet_name)
        df = pd.read_excel(file, sheet_name=sheet_name)
        col = str(style_df.data_df.columns[0])
        idx_start = [int(i) for i in style_df.index if idx_start_str in str(style_df.loc[i, col])][0]
        idx_end = [int(i) for i in style_df.index if idx_end_str in str(style_df.loc[i, col])][0]
        df.loc[idx_start:idx_end, 'indent'] = [int(i) for i in style_df[col].style.indent.loc[idx_start:idx_end]]

        df = df.dropna(subset=['indent'])
        tot = df.iloc[0]
        df['section'] = df.where(df.indent == 1).iloc[:, 0].ffill()
        df['vehicle_type'] = df.where(df.indent == 2).iloc[:, 0].ffill()

        if sheet_name == 'TrRoad_act':
            df['vehicle_subtype'] = df.where(df.indent == 3).iloc[:, 0]
            df = df.where(
                (df.indent == 3) | (df.vehicle_type == 'Powered 2-wheelers')
            ).dropna(how='all')
            df.loc[df.vehicle_type == 'Powered 2-wheelers', 'vehicle_subtype'] = 'Gasoline engine'
            df = (
                df
                .set_index(['section', 'vehicle_type', 'vehicle_subtype'])
                .drop([col, 'indent'], axis=1)
            )
        elif sheet_name == 'TrRoad_ene':
            df['vehicle_subtype'] = df.where(df.indent == 3).iloc[:, 0].ffill()
            # Deal with the annoying nature of the spreadsheets
            df['vehicle_type'] = df['vehicle_type'].str.split('(', expand=True)[0].str.strip()
            df['vehicle_subtype'] = df['vehicle_subtype'].str.split('(', expand=True)[0].str.strip()
            df.loc[df.vehicle_type == 'Powered 2-wheelers', 'vehicle_subtype'] = 'Gasoline engine'
            df['carrier'] = df.where(df.indent == 4).iloc[:, 0]
            df['carrier'] = df['carrier'].str.replace('of which ', '')
            df.loc[(df.vehicle_type == 'Powered 2-wheelers') & (df.indent == 2), 'carrier'] = 'petrol'
            df.loc[(df.vehicle_type == 'Powered 2-wheelers') & (df.indent == 3), 'carrier'] = 'biofuels'
            df.loc[(df.vehicle_subtype == 'Domestic') & (df.indent == 3), 'carrier'] = 'diesel'
            df.loc[(df.vehicle_subtype == 'International') & (df.indent == 3), 'carrier'] = 'diesel'
            df['carrier'] = df['carrier'].fillna(df.vehicle_subtype.replace({
                'Gasoline engine': 'petrol', 'Diesel oil engine': 'diesel',
                'Natural gas engine': 'natural_gas', 'LPG engine': 'lpg',
                'Battery electric vehicles': 'electricity', 'Plug-in hybrid electric': 'petrol'
            }))
            df = (
                df
                .where((df.indent > 2) | (df.vehicle_type == 'Powered 2-wheelers'))
                .dropna()
                .set_index(['section', 'vehicle_type', 'vehicle_subtype', 'carrier'])
                .drop([col, 'indent'], axis=1)
            )

            df = remove_of_which(df, 'diesel', 'biofuels')
            df = remove_of_which(df, 'petrol', 'biofuels')
            df = remove_of_which(df, 'petrol', 'electricity')
            df = remove_of_which(df, 'natural_gas', 'biogas')

        elif sheet_name == 'TrRail_ene' or sheet_name == 'TrRail_act':
            df['carrier'] = df.where(df.indent == 3).iloc[:, 0]
            df.loc[df.vehicle_type == 'Metro and tram, urban light rail', 'carrier'] = 'electricity'
            df.loc[df.vehicle_type == 'High speed passenger trains', 'carrier'] = 'electricity'
            carriers = {
                'Diesel oil (incl. biofuels)': 'diesel',
                'Electric': 'electricity',
                'Diesel oil': 'diesel'
            }
            df['carrier'] = df['carrier'].replace(carriers).fillna(df.vehicle_type.replace(carriers))
            df.loc[df.section == 'Freight transport', 'vehicle_type'] = 'Freight'

            df = (
                df
                .where((df.indent > 1) & (df.carrier.str.find('Conventional') == -1))
                .dropna()
                .set_index(['section', 'vehicle_type', 'carrier'])
                .drop([col, 'indent'], axis=1)
            )

        df = (
            df
            .assign(country_code=col.split(' - ')[0], unit=unit)
            .set_index(['country_code', 'unit'], append=True)
        )
        df.columns = df.columns.astype(int).rename('year')

        # After all this hardcoded cleanup, make sure numbers match up
        assert np.allclose(df.sum(), tot.loc[df.columns].astype(float))

        dfs.append(df)
    return pd.concat(dfs)


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
