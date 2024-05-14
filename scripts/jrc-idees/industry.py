import logging
from itertools import product
from multiprocessing import Pool
from typing import Callable, Literal, Union

import numpy as np
import pandas as pd
import xarray as xr
from eurocalliopelib import utils
from styleframe import StyleFrame

idx = pd.IndexSlice

SHEETS = ["ISI", "NFM", "CHI", "NMM", "PPA", "FBT", "TRE", "MAE", "TEL", "WWP", "OIS"]

LOGGER = logging.getLogger(__name__)

ENERGY_SHEET_COLORS = {
    "section": ["FFC00000"],
    "subsection": ["FF0070C0", "4f6228", "953735"],
    "end_use": ["984807"],
    "carrier_name": ["808080", "e46c0a", "000000", "dc9e9c", "d99694"],
    "skip_": [  # we have these to ensure we know what all sheet colours refer to
        "FF002060",  # general descriptor of sheet data, e.g. "Detailed split of energy consumption by subsector (ktoe)"
        "10253f",  # top row year references
    ],
}

ENERGY_SHEET_CARRIERS = {
    "electric|microwave|freeze|mechanical": "Electricity",
    "diesel": "Diesel oil (incl. biofuels)",
    "gas|thermal cooling|thermal connection": "Natural gas (incl. biogas)",
    "gas": "Natural gas (incl. biogas)",
}


def process_jrc_industry_data(
    paths_to_data: list[str],
    dataset: Literal["energy", "production"],
    threads: int,
    out_path: str,
):
    """Process human-readable JRC-IDEES Excel files into machine-readable datasets.
    This requires using the colour of Excel cell contents to classify rows.
    The classification of colours in `ENERGY_SHEET_COLORS` has been undertaken manually.

    Args:
        data_dir (str): Directory in which files are stored
        dataset (Literal[energy, production]): Data to process from those files.
        threads (int): Number of multi-processing threads to use.
        out_path (str): Path to which machine-readable data will be stored.
    """

    if dataset == "energy":
        processed_data = process_sheets(paths_to_data, threads, get_jrc_idees_energy)
        unit = "twh"
    elif dataset == "production":
        processed_data = process_sheets(paths_to_data, 1, get_jrc_idees_production)
        unit = "kt"

    processed_xr_data = df_to_xr(processed_data, unit)
    processed_xr_data.to_netcdf(out_path)


def df_to_xr(df: pd.DataFrame, unit: str) -> Union[xr.Dataset, xr.DataArray]:
    df.columns = df.columns.rename("year").astype(int)

    xr_data = df.stack().unstack("variable").to_xarray()

    country_code_mapping = utils.convert_valid_countries(
        xr_data.country_code.values, errors="ignore"
    )
    xr_data = utils.rename_and_groupby(
        xr_data, country_code_mapping, dim_name="country_code"
    )

    for var in xr_data.data_vars:
        xr_data[var] = xr_data[var].assign_attrs(units=unit)

    return xr_data


def process_sheets(
    data_filepaths: list[str], threads: int, processing_script: Callable
) -> pd.DataFrame:
    "Process energy sheet in data files across multiple threads"
    with Pool(int(threads)) as pool:
        dfs = pool.starmap(processing_script, product(SHEETS, data_filepaths))

    processed_df = pd.concat(dfs).sort_index()

    return processed_df


def get_jrc_idees_production(sheet_name: str, file: str) -> pd.DataFrame:
    xls = pd.ExcelFile(file)
    df = pd.read_excel(xls, sheet_name=sheet_name, index_col=0)
    start = df.filter(regex="Physical output", axis=0)
    end = df.filter(regex="Installed capacity", axis=0)
    df_processed = (
        df.loc[start.index[0] : end.index[0]]
        .iloc[1:-1]
        .dropna(how="all")
        .assign(
            country_code=df.index.name.split(":")[0],
            cat_name=df.index.name.split(": ")[1],
        )
        .rename_axis(index="produced_material")
    )
    df_processed.index = (
        df_processed.index.str.replace("(kt)", "", regex=False)
        .str.replace("(kt ", "(", regex=False)
        .str.strip()
    )
    df_processed = df_processed.assign(variable="production")
    return df_processed.set_index(["variable", "country_code", "cat_name"], append=True)


def get_jrc_idees_energy(sheet: str, file: str) -> pd.DataFrame:
    LOGGER.info(f"Processing file: {file}, sheet: {sheet}")
    xls = pd.ExcelFile(file)
    final_energy = _get_jrc_idees_energy_sheet(f"{sheet}_fec", xls)
    useful_energy = _get_jrc_idees_energy_sheet(f"{sheet}_ued", xls)
    return pd.concat(
        [final_energy, useful_energy], names=["variable"], keys=["final", "useful"]
    )


def _get_jrc_idees_energy_sheet(sheet_name: str, xls: str) -> pd.DataFrame:
    """
    This sheet needs to be parsed both based on the colour of the cell and the indent
    level of the text inside the cell.
    """
    style_df = StyleFrame.read_excel(xls, read_style=True, sheet_name=sheet_name)
    df = pd.read_excel(xls, sheet_name=sheet_name)
    column_names = str(style_df.data_df.columns[0])

    last_index_of_data = int(
        style_df[style_df[column_names].str.find("Market shares") > -1].item()
    )
    df = _assign_section_level_based_on_colour(
        style_df, df, column_names, last_index_of_data
    )
    df, total_to_check = _slice_on_indent(
        style_df, df, column_names, last_index_of_data
    )

    df = _rename_carriers(df)
    df = _assign_category_country_information(df, column_names)

    # Check that we haven't lost some data
    assert np.allclose(
        total_to_check.reindex(df.columns).astype(float), df.sum().astype(float)
    )

    return df.apply(utils.ktoe_to_twh)


def _assign_section_level_based_on_colour(
    style_df: StyleFrame, df: pd.DataFrame, column_names: str, last_index_of_data: int
) -> pd.DataFrame:
    all_colours = set(style_df[column_names].style.font_color.loc[:last_index_of_data])

    for section_level, colours in ENERGY_SHEET_COLORS.items():
        all_colours.difference_update(colours)
        if section_level == "skip_":
            continue
        idx = (
            style_df[column_names]
            .style.font_color.isin(colours)
            .loc[:last_index_of_data]
        )
        df.loc[idx[idx].index, section_level] = (
            style_df.loc[idx[idx].index, column_names]
            .astype(str)
            .where(lambda x: x != "nan")
        )
    if all_colours:
        raise AssertionError(
            f"Some rows in {df.columns[0]} have not been allocated to index levels as they have unexpected colour(s): {all_colours}"
        )
    return df


def _slice_on_indent(
    style_df: StyleFrame, df: pd.DataFrame, column_names: str, last_index_of_data: int
) -> tuple[pd.DataFrame, pd.Series]:
    df = df.loc[:last_index_of_data]
    df = df.assign(
        indent=style_df[column_names].style.indent.loc[:last_index_of_data].astype(int)
    )

    # indent of 1 tab == high-level summed data.
    # We keep this to check at the end that we haven't lost some data
    total_to_check = df[df.indent == 1].sum()

    df = df.dropna(
        subset=["section", "subsection", "end_use", "carrier_name"], how="all"
    )

    # All data between section/subsection names correspond to the preceding section/subsection name
    df["section"].ffill(inplace=True)
    df["subsection"].ffill(inplace=True)

    # When the indent decreases from one row to the next it signals that the data in that row
    # is an aggregation of a set of previous data.
    # We don't want these aggregated data rows in our final result.
    df = df.where(df.indent >= df.indent.shift(-1).fillna(df.indent)).dropna(how="all")

    return df, total_to_check


def _rename_carriers(df: pd.DataFrame) -> pd.DataFrame:
    # ASSUME: if end uses / subsections without end use disaggregation do not have carrier disaggregation,
    # then their name gives us enough of a clue as to the energy carrier their energy values refer to.
    # E.g., microwave == electricity, thermal == gas powered, mention of "diesel" == diesel oil.
    for carrier_search_string, carrier_group in ENERGY_SHEET_CARRIERS.items():
        df.loc[
            df.end_use.fillna(df.subsection)
            .str.lower()
            .str.contains(carrier_search_string, regex=True, na=False),
            "carrier_name",
        ] = carrier_group
    # ASSUME: if we cannot infer anything from the end use / subsection name, the carrier is electricity
    to_fill = df[df.carrier_name.isnull()]
    if not to_fill.empty:
        LOGGER.info(
            f"Filling the following rows in {df.columns[0]} with the carrier 'Electricity': {set(to_fill.iloc[:, 0].values)} "
        )
        df.carrier_name.fillna("Electricity", inplace=True)
    return df


def _assign_category_country_information(
    df: pd.DataFrame, column_names: str
) -> pd.DataFrame:
    index = ["section", "subsection", "carrier_name", "country_code", "cat_name"]
    return (
        df.assign(
            cat_name=column_names.split(": ")[1].split(" / ")[0],
            country_code=column_names.split(": ")[0],
        )
        .set_index(index)
        .drop([column_names, "indent", "end_use"], axis="columns")
        .sum(level=index)
    )


if __name__ == "__main__":
    process_jrc_industry_data(
        paths_to_data=snakemake.input.data,
        dataset=snakemake.wildcards.dataset,
        threads=snakemake.threads,
        out_path=snakemake.output[0],
    )
