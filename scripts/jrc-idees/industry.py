import glob
import os

import numpy as np
import pandas as pd
from styleframe import StyleFrame

from eurocalliopelib import utils

idx = pd.IndexSlice

SHEETS = ['ISI', 'NFM', 'CHI', 'NMM', 'PPA', 'FBT', 'TRE', 'MAE', 'TEL', 'WWP', 'OIS']


def process_jrc_industry_data(data_dir, dataset, out_path):
    data_filepaths = glob.glob(os.path.join(data_dir, "*.xlsx"))
    if dataset == "energy":
        processed_data = pd.concat(
            [
                pd.concat(
                    get_jrc_idees_energy_sheet(f'{sheet}_{flow}', data_filepaths)
                    for sheet in SHEETS
                )
                for flow in ["fec", "ued"]
            ],
            names=['energy'], keys=['consumption', 'demand']
        )
        processed_data = processed_data.apply(utils.ktoe_to_twh)
        processed_data.index = processed_data.index.set_levels(['twh'], level='unit')
    elif dataset == "production":
        processed_data = pd.concat(get_jrc_idees_production(sheet, data_filepaths) for sheet in SHEETS)
    processed_data.stack('year').to_csv(out_path)


def get_jrc_idees_energy_sheet(sheet_name, files):
    dfs = []
    for file in files:
        style_df = StyleFrame.read_excel(file, read_style=True, sheet_name=sheet_name)
        df = pd.read_excel(file, sheet_name=sheet_name)
        col = str(style_df.data_df.columns[0])
        idx_end = [int(i) for i in style_df.index if 'Market shares' in str(style_df.loc[i, col])][0]
        colors = {
            'section': ['FFC00000'],
            'subsection': ['FF0070C0', '4f6228', '953735'],
            'end_use': ['984807'],
            'carrier_name': ['808080', 'e46c0a', '000000', 'dc9e9c']
        }

        for k, v in colors.items():
            idx = style_df[col].style.font_color.isin(v).loc[:idx_end]
            df.loc[idx[idx == True].index, k] = [
                str(i) if str(i) != 'nan' else np.nan
                for i in style_df.loc[idx[idx == True].index, col]
            ]
        df.loc[:idx_end, 'indent'] = [int(i) for i in style_df[col].style.indent.loc[:idx_end]]

        # To check at the end that we haven't lost some data
        tot = df[df.indent == 1].sum()

        df = df.dropna(subset=['section', 'subsection', 'end_use', 'carrier_name'], how='all')
        df['section'] = df.section.ffill()
        df['subsection'] = df.subsection.ffill()
        df.loc[
            df.end_use.str.lower().str.contains('electric|microwave|freeze|mechanical', regex=True, na=False),
            'carrier_name'
        ] = 'Electricity'
        df.loc[
            df.subsection.str.lower().str.contains('diesel', regex=True, na=False),
            'carrier_name'
        ] = 'Diesel oil (incl. biofuels)'
        df.loc[
            df.end_use.str.lower().str.contains('gas|thermal cooling|thermal connection', regex=True, na=False),
            'carrier_name'
        ] = 'Natural gas (incl. biogas)'
        df.loc[
            df.end_use.str.lower().str.contains('gas', regex=True, na=False),
            'carrier_name'
        ] = 'Natural gas (incl. biogas)'
        df.loc[df.end_use.isnull(), 'carrier_name'] = (
            df.loc[df.end_use.isnull(), 'carrier_name'].fillna('Electricity')
        )

        df = df[df.indent >= df.indent.shift(-1).fillna(df.indent)]

        df = (
            df
            .assign(cat_name=col.split(': ')[1].split(' / ')[0], country_code=col.split(': ')[0], unit='ktoe')
            .set_index(['section', 'subsection', 'carrier_name', 'country_code', 'cat_name', 'unit'])
            .drop([col, 'indent', 'end_use'], axis=1)
            .groupby(level=[0, 1, 2, 3, 4, 5]).sum()
        )
        # Check that we haven't lost some data
        assert np.allclose(tot.reindex(df.columns).astype(float), df.sum().astype(float))

        dfs.append(df)

    return pd.concat(dfs)


def get_jrc_idees_production(sheet_name, files):
    dfs = []
    for file in files:
        df = pd.read_excel(file, sheet_name=sheet_name, index_col=0)
        start = df.filter(regex='Physical output', axis=0)
        end = df.filter(regex='Installed capacity', axis=0)

        df_prod = (
            df
            .loc[start.index[0]:end.index[0]]
            .iloc[1:-1]
            .dropna(how='all')
            .assign(country_code=df.index.name.split(':')[0], cat_name=df.index.name.split(': ')[1], unit='kt')
            .rename_axis(index='produced_material')
            .set_index(['country_code', 'cat_name', 'unit'], append=True)
        )

        dfs.append(df_prod)

    return pd.concat(dfs)


if __name__ == "__main__":
    process_jrc_industry_data(
        data_dir=snakemake.input.unprocessed_data,
        dataset=snakemake.wildcards.dataset,
        out_path=snakemake.output[0]
    )
