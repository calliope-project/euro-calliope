from pathlib import Path

import numpy as np
import pandas as pd

from eurocalliopelib import utils

idx = pd.IndexSlice

END_USES = {
    'Space heating': 'space_heating',
    'Space cooling': 'cooling',
    'Hot water': 'water_heating',
    'Catering': 'cooking'
}
CARRIER_NAMES = {
    'Advanced electric heating': 'electricity',
    'Biomass and wastes': 'biofuel',
    'Conventional electric heating': 'electricity',
    'Conventional gas heaters': 'gas',
    'Derived heat': 'heat',
    'Electric space cooling': 'electricity',
    'Electricity': 'electricity',
    'Electricity in circulation and other use': 'electricity',
    'Gas heat pumps': 'gas',
    'Gas/Diesel oil incl. biofuels (GDO)': 'oil',
    'Gases incl. biogas': 'gas',
    'Geothermal energy': 'renewable_heat',
    'Liquified petroleum gas (LPG)': 'oil',
    'Solar': 'renewable_heat',
    'Solids': 'solid_fossil'
}


def process_jrc_tertiary_data(data_dir, out_path):
    data_filepaths = list(Path(data_dir).glob("*.xlsx"))
    processed_data = pd.concat([get_tertiary_sector_data(file) for file in data_filepaths])
    processed_data = processed_data.apply(utils.ktoe_to_twh)

    processed_da = processed_data.rename("jrc-idees-tertiary-twh").to_xarray()

    country_code_mapping = utils.convert_valid_countries(processed_da.country_code.values)
    processed_da = utils.rename_and_groupby(processed_da, country_code_mapping, dim="country_code")

    processed_da.assign_attrs(unit="twh").to_netcdf(out_path)


def get_tertiary_sector_data(file):
    xls = pd.ExcelFile(file)
    df_consumption = pd.read_excel(xls, sheet_name='SER_hh_fec', index_col=0)
    df_demand = pd.read_excel(xls, sheet_name='SER_hh_tes', index_col=0)
    df_summary = pd.read_excel(xls, sheet_name='SER_summary', index_col=0)

    df_consumption = clean_df(df_consumption, 'consumption')
    df_demand = clean_df(df_demand, 'demand')

    df = pd.concat([df_consumption, df_demand]).sort_index()

    df = add_electricity_use(df, df_summary)
    assert np.allclose(
        df.xs('consumption', level='energy').sum(),
        df_summary.loc['Energy consumption by fuel - Eurostat structure (ktoe)'].astype(float)
    )

    return df.stack()


def clean_df(df, energy_type):
    # The index is a flattened multi-level index, where the end-use is the highest aggregation level
    # in the hierarchy. We identify this point in the hierarchy and use ffill to match lower aggregation levels
    # with the top-level end-use.
    df.loc[df.index.isin(END_USES.keys()), 'end_use'] = list(END_USES.keys())
    df.end_use = df.end_use.ffill()

    df = (
        df
        .dropna()
        .set_index('end_use', append=True)
        .drop(END_USES.keys(), level=0)  # we remove the top-level end-use aggregation level
        .groupby([CARRIER_NAMES, END_USES], level=[0, 1]).sum()
        .assign(
            country_code=df.index.names[0].split(' - ')[0],
            energy=energy_type
        )
        .set_index(['country_code', 'energy'], append=True)
        .rename_axis(columns='year', index=['carrier_name', 'end_use', 'country_code', 'energy'])
    )
    return df


def add_electricity_use(df, df_summary):
    """End-use electricity consumption is added equally to both 'demand' and 'consumption'"""
    df_elec = (
        df_summary
        .loc['Energy consumption by end-uses (ktoe)':'Shares of energy consumption in end-uses (in %)']
        .loc['Specific electricity uses']
        .rename_axis(index='year')
    )
    new_idx = pd.MultiIndex.from_product(
        [["electricity"], ["end_use_electricity"], df.index.levels[2], ["consumption", "demand"]],
        names=df.index.names
    )
    df = df.append(pd.concat([df_elec, df_elec], keys=new_idx, axis=1).T).sort_index()
    return df


if __name__ == "__main__":
    process_jrc_tertiary_data(
        data_dir=snakemake.input.unprocessed_data,
        out_path=snakemake.output[0]
    )
