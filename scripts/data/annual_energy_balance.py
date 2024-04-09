from enum import Enum
from string import digits

import pandas as pd
from eurocalliopelib import utils

idx = pd.IndexSlice


class CAT_CODE(Enum):
    FINAL_CONSUMPTION_HOUSEHOLD_CATEGORY = "FC_OTH_HH_E"
    FINAL_CONSUMPTION_INDUSTRY_CATEGORY = "FC_IND_E"
    FINAL_CONSUMPTION_OTHER_SECTORS_COMMERCIAL_PUBLIC_SERVICES = "FC_OTH_CP_E"


def generate_annual_energy_balance_nc(
    path_to_energy_balance: str,
    path_to_cat_names: str,
    path_to_carrier_names: str,
    path_to_ch_excel: str,
    path_to_ch_industry_excel: str,
    path_to_result: str,
    first_year: int,
) -> None:
    """
    Open a TSV file and reprocess it into a xarray dataset, including long names for
    Eurostat codes.
    Switzerland is not included in Eurostat, so we splice in data from their govt.
    statistics.
    """
    # Names for each consumption category/sub-category and carriers have been prepared by hand
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)

    df = pd.read_csv(
        path_to_energy_balance,
        delimiter="\t",
        index_col=0,
        na_values=[":", ": ", ": z"],
    )
    df.index = (
        df.index.str.split(",", expand=True).rename([
            "cat_code",
            "carrier_code",
            "unit",
            "country",
        ])  # comes as 'nrg_bal,siec,unit,geo\\time'
    )
    not_countries = [c for c in df.reset_index().country.unique() if len(c) > 2] + [
        "XK"
    ]
    df = (
        df.drop(axis=0, level="country", labels=not_countries)
        .reset_index(level="country")
        .assign(country=lambda df: df.country.map(utils.convert_country_code))
        .set_index("country", append=True)
    )
    df.columns = df.columns.astype(int).rename("year")
    df = df.loc[idx[cat_names.index, carrier_names.index, "TJ"], :].dropna(how="all")
    df = df.sort_index(axis=1).loc[:, first_year:]

    tdf = df.stack()

    # Add CH energy use (only covers a subset of sectors and carriers, but should be enough)
    ch_energy_use_tdf = add_ch_energy_balance(
        path_to_ch_excel, path_to_ch_industry_excel, index_levels=tdf.index.names
    )
    tdf = pd.concat([tdf, ch_energy_use_tdf]).sort_index(axis=0)

    # TODO treat missing values if necessary

    tdf.rename("value").to_csv(path_to_result)


def add_ch_energy_balance(path_to_ch_excel, path_to_ch_industry_excel, index_levels):
    household_sheet = "T17a"
    industry_sheet = "T17b"
    other_sectors_sheet = "T17c"

    ch_hh_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel,
        household_sheet,
        skipfooter=9,
        cat_code=CAT_CODE.FINAL_CONSUMPTION_HOUSEHOLD_CATEGORY,
    )
    ch_ind_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel,
        industry_sheet,
        skipfooter=12,
        cat_code=CAT_CODE.FINAL_CONSUMPTION_INDUSTRY_CATEGORY,
    )
    ch_ser_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel,
        other_sectors_sheet,
        skipfooter=12,
        cat_code=CAT_CODE.FINAL_CONSUMPTION_OTHER_SECTORS_COMMERCIAL_PUBLIC_SERVICES,
    )

    ch_waste_energy_use = get_ch_waste_consumption(path_to_ch_excel)
    ch_industry_subsector_energy_use = get_ch_industry_energy_balance(
        path_to_ch_industry_excel
    )
    ch_transport_energy_use = get_ch_transport_energy_balance(path_to_ch_excel)

    ch_energy_use_tdf = pd.concat([
        df.reset_index("year")
        .assign(country="CHE", unit="TJ")
        .set_index(["year", "country", "unit"], append=True)
        .squeeze()
        .reorder_levels(index_levels)
        for df in [
            ch_hh_energy_use,
            ch_ind_energy_use,
            ch_ser_energy_use,
            ch_waste_energy_use,
            ch_industry_subsector_energy_use,
            ch_transport_energy_use,
        ]
    ])

    return ch_energy_use_tdf


def get_ch_energy_balance_sheet(path_to_excel, sheet, skipfooter, cat_code):
    ch_energy_carriers = {
        "Erdölprodukte": "O4000XBIO",
        "Elektrizität": "E7000",
        "Gas": "G3000",
        "Kohle": "C0000X0350-0370",
        "Holzenergie": "R5110-5150_W6000RI",
        "Fernwärme": "H8000",
        "Industrieabfälle": "W6100_6220",
        "Übrige erneuerbare Energien": "RA000",
        "Total\n= %": "TOTAL",
    }
    # Footnote labels lead to some strings randomly ending in numbers; we remove them here
    remove_digits = str.maketrans("", "", digits)
    df = (
        pd.read_excel(
            path_to_excel,
            skiprows=6,
            skipfooter=skipfooter,
            index_col=0,
            sheet_name=sheet,
            dtype="float",
            na_values=["-"],
            header=[0, 1, 2, 3, 4],
        )
        .xs("TJ", level=-1, axis=1)  # Ignore columns giving % use
        .drop(
            columns=["Erdölprodukte", "Erdölprodukte1"], level=0, errors="ignore"
        )  # ignore the column giving subset of oil use which is light oil
        .rename_axis(index="year")
    )
    df.columns = (
        df.columns.get_level_values(0)
        .str.translate(remove_digits)
        .map(ch_energy_carriers)
        .rename("carrier_code")
    )

    return df.assign(cat_code=cat_code).set_index("cat_code", append=True).stack()


def get_ch_waste_consumption(path_to_excel):
    """
    In a different sheet in the CH GEST dataset, get data on the consumed quantity of
    waste burned in WtE plants, ignoring the small (~2-3%) quantity of fossil fuels
    also consumed in WtE plants to kickstart the process.
    ASSUME: Small quantity (~2-3%) of fossil fuels consumed in Swiss WtE plants can be ignored.
    """
    category_code = "TI_EHG_E"
    carrier_code = "W6100_6220"
    sheet_name = "T27"

    waste_stream_gwh = pd.read_excel(
        path_to_excel,
        sheet_name=sheet_name,
        skiprows=5,
        index_col=0,
        header=[0, 1],
        skipfooter=8,
    )[("Consommation d'énergie (GWh)", "Ordures")]
    waste_stream_tj = waste_stream_gwh.apply(utils.gwh_to_tj)
    waste_stream_tdf = (
        waste_stream_tj.to_frame(carrier_code)  # carrier code
        .rename_axis(index="year", columns="carrier_code")
        .assign(cat_code=category_code)  # cat code
        .set_index("cat_code", append=True)
        .stack()
    )
    return waste_stream_tdf


def get_ch_transport_energy_balance(path_to_excel):
    carriers = {
        "Gas übriger Vekehr": ("G3000", "FC_TRA_ROAD_E"),
        "Übrige erneuerbare Energien": ("R5220B", "FC_TRA_ROAD_E"),
        "davon Benzin": ("O4652XR5210B", "FC_TRA_ROAD_E"),
        "davon Diesel": ("O4671XR5220B", "FC_TRA_ROAD_E"),
        "davon Flugtreibstoffe": ("O4000XBIO", "INTAVI"),
        "davon Bahnen": ("E7000", "FC_TRA_RAIL_E"),
        "davon Strasse": ("E7000", "FC_TRA_ROAD_E"),
    }
    # ASSUME "davon Non-Road" is not included in electrified transport

    df = pd.read_excel(
        path_to_excel,
        skiprows=6,
        skipfooter=12,
        index_col=0,
        sheet_name="T17e",
        dtype="float",
        na_values=["-"],
        header=[0, 1, 2, 3, 4],
    ).xs("TJ", level=-1, axis=1)

    # Footnote labels lead to some strings randomly ending in numbers; we remove them here
    remove_digits = str.maketrans("", "", digits)
    # carrier names span across two column levels, which we merge with fillna

    def carrier_name_func(index):
        return (
            df.columns.to_frame()
            .iloc[:, index]
            .str.translate(remove_digits)
            .map(carriers)
        )

    df.columns = carrier_name_func(0).fillna(carrier_name_func(1)).values

    df = (
        df
        .groupby(axis=1, level=0)
        .sum()
        .rename_axis(index="year", columns="carrier_code")
        .T
    )
    df.index = pd.MultiIndex.from_tuples(df.index, names=('carrier_code', 'cat_code'))
    return df.stack()


def get_ch_industry_energy_balance(path_to_excel):
    ch_subsectors = {
        "1 Nahrg.": "FC_IND_FBT_E",  # 'Food, beverages & tobacco',
        "2 Textil": "FC_IND_TL_E",  # 'Textile & leather',
        "3 Papier": "FC_IND_PPP_E",  # 'Paper, pulp & printing',
        "4 Chemie": "FC_IND_CPC_E",  # 'Chemical & petrochemical',
        "5 Zement": "FC_IND_NMM_E",  # 'Non-metallic minerals',
        "6 andere": "FC_IND_NMM_E",  # 'Non-metallic minerals',
        "7 Metall": "FC_IND_IS_E",  # 'Iron & steel',
        "8 NE": "FC_IND_NFM_E",  # 'Non-ferrous metals',
        "9 Metall": "FC_IND_MAC_E",  # 'Machinery',
        "10 Masch": "FC_IND_MAC_E",  # 'Machinery',
        "11 and.": "FC_IND_NSP_E",
        # 'Not elsewhere specified (industry)', 'Wood & wood products, Mining & quarrying, Transport equipment
        "12 Bau": "FC_IND_CON_E",  # 'Construction'
    }

    ch_carriers = {  # first row in which carriers are defined in the file
        25: "E7000",  # 'electricity',
        52: "O4000XBIO",  # 'oil',
        79: "G3000",  # 'gas',
        105: "C0000X0350-0370",  # 'solid_fuel',
        128: "W6100_6220",  # 'waste',
        151: "O4000XBIO",  # 'oil',
        191: "H8000",  # 'heat',  # purchased
        228: "R5110-5150_W6000RI",  # 'biofuel'
    }

    column_names = (
        pd.read_excel(
            path_to_excel,
            sheet_name="Überblick_tot",
            skiprows=5,
            usecols="A,E:P",
            nrows=0,
            header=0,
        )
        .rename(columns=ch_subsectors)
        .iloc[:, 1:]  # ignore index column
    )
    column_names = ["year"] + list(column_names)

    return (
        pd.concat(
            [
                read_industry_subsector(
                    path_to_excel, first_row, column_names, carrier_name
                )
                for first_row, carrier_name in ch_carriers.items()
            ],
            axis=0,
        )
        .groupby(axis=0, level=[0, 1])  # group carriers
        .sum()
        .stack()
        .rename("value")
    )


def read_industry_subsector(
    path: str, first_row: int, column_names: list[str], carrier_code: str
) -> pd.DataFrame:
    return (
        pd.read_excel(
            path,
            sheet_name="Überblick_tot",
            skiprows=first_row - 1,
            usecols="A,E:P",
            nrows=10,
            header=None,
            names=column_names,
            index_col=0,
        )
        .rename_axis(columns="cat_code")
        .rename(index=lambda year: str(year)[:4])
        .rename(
            columns=lambda sector: sector.split(".")[0]
        )  # pandas adds dots and numbers to duplicates, remove
        .assign(carrier_code=carrier_code)
        .reset_index()
        .set_index(["carrier_code", "year"])
        .groupby(axis=1, level=0)  # group sectors
        .sum()
    )


if __name__ == "__main__":
    generate_annual_energy_balance_nc(
        path_to_energy_balance=snakemake.input.energy_balance,
        path_to_ch_excel=snakemake.input.ch_energy_balance,
        path_to_ch_industry_excel=snakemake.input.ch_industry_energy_balance,
        path_to_cat_names=snakemake.input.cat_names,
        path_to_carrier_names=snakemake.input.carrier_names,
        first_year=snakemake.params.first_year,
        path_to_result=snakemake.output[0],
    )
