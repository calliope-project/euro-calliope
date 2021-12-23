from pathlib import Path
from multiprocessing import Pool
from itertools import product

import numpy as np
import pandas as pd
from styleframe import StyleFrame

from eurocalliopelib import utils

idx = pd.IndexSlice

SHEETS = ['ISI', 'NFM', 'CHI', 'NMM', 'PPA', 'FBT', 'TRE', 'MAE', 'TEL', 'WWP', 'OIS']

ENERGY_SHEET_COLORS = {
    'section': ['FFC00000'],
    'subsection': ['FF0070C0', '4f6228', '953735'],
    'end_use': ['984807'],
    'carrier_name': ['808080', 'e46c0a', '000000', 'dc9e9c']
}

ENERGY_SHEET_CARRIERS = {
    "electric|microwave|freeze|mechanical": "Electricity",
    "diesel": "Diesel oil (incl. biofuels)",
    "gas|thermal cooling|thermal connection": "Natural gas (incl. biogas)",
    "gas": "Natural gas (incl. biogas)",
}


def process_jrc_industry_data(data_dir, dataset, threads, out_path):
    data_filepaths = list(Path(data_dir).glob("*.xlsx"))
    if dataset == "energy":
        processed_data = process_energy(data_filepaths, threads)
    elif dataset == "production":
        processed_data = process_production(data_filepaths)

    processed_data.columns = processed_data.columns.rename("year").astype(int)
    processed_data.stack().to_csv(out_path)


def process_energy(data_filepaths, threads):
    with Pool(int(threads)) as pool:
        demand_consumption = pool.starmap(
            get_jrc_idees_demand_consumption,
            product(SHEETS, data_filepaths)
        )

    processed_data = pd.concat(demand_consumption).sort_index()
    processed_data = processed_data.apply(utils.ktoe_to_twh)
    processed_data.index = processed_data.index.set_levels(['twh'], level='unit')

    return processed_data


def process_production(data_filepaths):
    all_dfs = []
    for file in data_filepaths:
        xls = pd.ExcelFile(file)
        all_dfs.extend([get_jrc_idees_production_sheet(sheet, xls) for sheet in SHEETS])

    return pd.concat(all_dfs)


def get_jrc_idees_demand_consumption(sheet, file):
    print(sheet, file)
    xls = pd.ExcelFile(file)
    consumption = get_jrc_idees_energy_sheet(f"{sheet}_fec", xls)
    demand = get_jrc_idees_energy_sheet(f"{sheet}_ued", xls)
    return pd.concat(
        [consumption, demand],
        names=['energy'],
        keys=['consumption', 'demand']
    )


def get_jrc_idees_energy_sheet(sheet_name, xls):
    """
    This sheet needs to be parsed both based on the colour of the cell and the indent
    level of the text inside the cell.
    """
    style_df = StyleFrame.read_excel(xls, read_style=True, sheet_name=sheet_name)
    df = pd.read_excel(xls, sheet_name=sheet_name)
    column_names = str(style_df.data_df.columns[0])

    last_index_of_data = int(style_df[style_df[column_names].str.find('Market shares') > -1].item())
    df = assign_section_level_based_on_colour(style_df, df, column_names, last_index_of_data)
    df, total_to_check = slice_on_indent(style_df, df, column_names, last_index_of_data)
    df = rename_carriers(df)
    df = assign_category_country_unit_information(df, column_names)

    # Check that we haven't lost some data
    assert np.allclose(
        total_to_check.reindex(df.columns).astype(float),
        df.sum().astype(float)
    )

    return df


def assign_section_level_based_on_colour(style_df, df, column_names, last_index_of_data):

    for section_level, colours in ENERGY_SHEET_COLORS.items():
        idx = style_df[column_names].style.font_color.isin(colours).loc[:last_index_of_data]
        df.loc[idx[idx].index, section_level] = (
            style_df.loc[idx[idx].index, column_names].astype(str).where(lambda x: x != 'nan')
        )
    return df


def slice_on_indent(style_df, df, column_names, last_index_of_data):
    df = df.loc[:last_index_of_data]
    df = df.assign(
        indent=style_df[column_names].style.indent.loc[:last_index_of_data].astype(int)
    )

    # indent of 1 tab == high-level summed data.
    # We keep this to check at the end that we haven't lost some data
    total_to_check = df[df.indent == 1].sum()

    df = df.dropna(subset=['section', 'subsection', 'end_use', 'carrier_name'], how='all')

    # All data between section/subsection names correspond to the preceding section/subsection name
    df["section"].ffill(inplace=True)
    df["subsection"].ffill(inplace=True)

    # When the indent decreases from one row to the next it signals that the data in that row
    # is an aggregation of a set of previous data.
    # We don't want these aggregated data rows in our final result.
    df = df.where(df.indent >= df.indent.shift(-1).fillna(df.indent)).dropna(how="all")

    return df, total_to_check


def rename_carriers(df):
    for carrier_search_string, carrier_group in ENERGY_SHEET_CARRIERS.items():
        df.loc[
            df.end_use.str.lower().str.contains(
                carrier_search_string, regex=True, na=False
            ),
            'carrier_name'
        ] = carrier_group
    df.loc[df.end_use.isnull(), 'carrier_name'] = (
        df.loc[df.end_use.isnull(), 'carrier_name'].fillna('Electricity')
    )
    return df


def assign_category_country_unit_information(df, column_names):
    index = ['section', 'subsection', 'carrier_name', 'country_code', 'cat_name', 'unit']
    return (
        df
        .assign(
            cat_name=column_names.split(': ')[1].split(' / ')[0],
            country_code=column_names.split(': ')[0],
            unit='ktoe'
        )
        .set_index(index)
        .drop([column_names, 'indent', 'end_use'], axis="columns")
        .sum(level=index)
    )


def get_jrc_idees_production_sheet(sheet_name, xls):
    df = pd.read_excel(xls, sheet_name=sheet_name, index_col=0)
    start = df.filter(regex='Physical output', axis=0)
    end = df.filter(regex='Installed capacity', axis=0)

    return (
        df
        .loc[start.index[0]:end.index[0]]
        .iloc[1:-1]
        .dropna(how='all')
        .assign(
            country_code=df.index.name.split(':')[0],
            cat_name=df.index.name.split(': ')[1],
            unit='kt'
        )
        .rename_axis(index='produced_material')
        .set_index(['country_code', 'cat_name', 'unit'], append=True)
    )


if __name__ == "__main__":
    process_jrc_industry_data(
        data_dir=snakemake.input.unprocessed_data,
        dataset=snakemake.wildcards.dataset,
        threads=snakemake.threads,
        out_path=snakemake.output[0]
    )
