
import pandas as pd

from eurocalliopelib import utils

CH_ADMINISTRATIVE_UNIT_TO_NUTS3 = {
    'AG': 'CH033', 'AI': 'CH054', 'AR': 'CH053', 'BE': 'CH021',
    'BL': 'CH032', 'BS': 'CH031', 'FR': 'CH022', 'GE': 'CH013',
    'GL': 'CH051', 'GR': 'CH056', 'JU': 'CH025', 'LU': 'CH061',
    'NE': 'CH024', 'NW': 'CH065', 'OW': 'CH064', 'SG': 'CH055',
    'SH': 'CH052', 'SO': 'CH023', 'SZ': 'CH063', 'TG': 'CH057',
    'TI': 'CH070', 'UR': 'CH062', 'VD': 'CH011', 'VS': 'CH012',
    'ZG': 'CH066', 'ZH': 'CH040'
}


def ch_gva(path_to_ch_gva_excel: str, path_to_output: str):
    administrative_units = pd.ExcelFile(path_to_ch_gva_excel).sheet_names
    ch_gva = pd.concat(
        [extract_gva_from_excel(path_to_ch_gva_excel, unit) for unit in administrative_units],
        names=["id", "year"], keys=administrative_units
    ).unstack("year")

    ch_gva.index = ch_gva.index.map(CH_ADMINISTRATIVE_UNIT_TO_NUTS3)
    ch_gva.columns = utils.to_numeric(ch_gva.columns).astype(int).rename("year")

    ch_gva.stack().to_csv(path_to_output)


def extract_gva_from_excel(path_to_ch_gva_excel, administrative_unit):
    ch_gva = pd.read_excel(
        path_to_ch_gva_excel, sheet_name=administrative_unit, index_col=0, skiprows=3, skipfooter=16,
        usecols="A,C:L"
    )
    # ASSUME: the subsectors 'G-J', 'K-N', and 'O-U' constitute commercial subsectors.
    # This matches the assumption made for equivalent eurostat data.
    return ch_gva.loc[["GHIJ", "K", "LMNRS", "O", "T"]].sum()


if __name__ == "__main__":
    ch_gva(
        path_to_ch_gva_excel=snakemake.input.ch_gva_excel,
        path_to_output=snakemake.output[0]
    )
