from string import digits

import pandas as pd
import numpy as np
import xarray as xr

from eurocalliopelib import utils

idx = pd.IndexSlice


def generate_annual_energy_balance_nc(
    path_to_input,
    path_to_cat_names,
    path_to_carrier_names,
    path_to_ch_excel,
    path_to_ch_industry_excel,
    path_to_result,
    countries,
):
    """
    Open a TSV file and reprocess it into a xarray dataset, including long names for
    Eurostat codes.
    Switzerland is not included in Eurostat, so we splice in data from their govt.
    statistics.
    """
    # Names for each consumption category/sub-category and carriers have been prepared by hand
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)
    country_codes = [utils.convert_country_code(i, output="alpha2_eurostat") for i in countries]

    df = pd.read_csv(path_to_input, delimiter="\t", index_col=0)
    df.index = df.index.str.split(",", expand=True).rename(
        ["cat_code", "carrier_code", "unit", "country"]
    )  # comes as 'nrg_bal,siec,unit,geo\\time'
    df.columns = df.columns.astype(int).rename("year")
    df = df.transform(utils.to_numeric)
    df = df.reorder_levels(["cat_code", "carrier_code", "unit", "country"])
    df = df.loc[
        idx[cat_names.index, carrier_names.index, "TJ", country_codes], :
    ].dropna(how="all")

    tdf = df.stack()

    # Add CH energy use (only covers a subset of sectors and carriers, but should be enough)
    ch_energy_use_tdf = add_ch_energy_balance(
        path_to_ch_excel, path_to_ch_industry_excel, index_levels=tdf.index.names
    )
    tdf = pd.concat([tdf, ch_energy_use_tdf]).sort_index(axis=0)

    tdf.to_csv(path_to_result)


def add_ch_energy_balance(path_to_ch_excel, path_to_ch_industry_excel, index_levels):
    ch_hh_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, "T17a", skipfooter=9, cat_code="FC_OTH_HH_E"
    )
    ch_ind_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, "T17b", skipfooter=12, cat_code="FC_IND_E"
    )
    ch_ser_energy_use = get_ch_energy_balance_sheet(
        path_to_ch_excel, "T17c", skipfooter=12, cat_code="FC_OTH_CP_E"
    )

    ch_waste_energy_use = get_ch_waste_consumption(path_to_ch_excel)
    ch_industry_subsector_energy_use = get_ch_industry_energy_balance(
        path_to_ch_industry_excel
    )
    ch_transport_energy_use = get_ch_transport_energy_balance(path_to_ch_excel)

    ch_energy_use_tdf = pd.concat(
        [
            i.reset_index("year")
            .assign(country="CH", unit="TJ")
            .set_index(["year", "country", "unit"], append=True)
            .reorder_levels(index_levels)
            for i in [
                ch_hh_energy_use,
                ch_ind_energy_use,
                ch_ser_energy_use,
                ch_waste_energy_use,
                ch_industry_subsector_energy_use,
                ch_transport_energy_use,
            ]
        ]
    )

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
    _df = (
        pd.read_excel(
            path_to_excel,
            skiprows=6,
            skipfooter=skipfooter,
            index_col=0,
            sheet_name=sheet,
            header=[0, 1, 2, 3, 4],
        )
        # Ignore columns giving % use
        .xs("TJ", level=-1, axis=1)
        # ignore the column giving subset of oil use which is light oil
        .iloc[:, [i for i in range(10) if i != 1]]
    )
    _df.columns = (
        _df.columns.get_level_values(0)
        .str.translate(remove_digits)
        .map(ch_energy_carriers)
        .rename("carrier_code")
    )
    _df.index.rename("year", inplace=True)

    _df = (
        _df.apply(utils.to_numeric)
        .astype(float)
        .assign(cat_code=cat_code)
        .set_index("cat_code", append=True)
        .stack()
    )

    return _df


def get_ch_waste_consumption(path_to_excel):
    """
    In a different sheet in the CH GEST dataset, get data on the consumed quantity of
    waste burned in WtE plants, ignoring the small (~2-3%) quantity of fossil fuels
    also consumed in WtE plants to kickstart the process.
    """
    waste_stream_gwh = pd.read_excel(
        path_to_excel,
        sheet_name="T27",
        skiprows=5,
        index_col=0,
        header=[0, 1],
        skipfooter=8,
    )[("Consommation d'énergie (GWh)", "Ordures")]
    waste_stream_tj = waste_stream_gwh.apply(utils.gwh_to_tj)
    waste_stream_tdf = (
        waste_stream_tj.to_frame("W6100_6220")  # carrier code
        .rename_axis(index="year", columns="carrier_code")
        .assign(cat_code="TI_EHG_E")  # cat code
        .set_index("cat_code", append=True)
        .stack()
    )
    return waste_stream_tdf


def get_ch_transport_energy_balance(path_to_excel):
    carriers = {
        "Elektrizität": "E7000",
        "Gas übriger Vekehr": "G3000",
        "Übrige erneuerbare Energien": "R5220B",
        "davon Benzin": "O4652XR5210B",
        "davon Diesel": "O4671XR5220B",
        "davon Flugtreibstoffe": "O4000XBIO",
    }
    categories = {
        "O4652XR5210B": "FC_TRA_ROAD_E",
        "O4671XR5220B": "FC_TRA_ROAD_E",
        "R5220B": "FC_TRA_ROAD_E",
        "G3000": "FC_TRA_ROAD_E",
        "E7000": "FC_TRA_RAIL_E",
        "O4000XBIO": "INTAVI",
    }
    _df = pd.read_excel(
        path_to_excel,
        skiprows=6,
        skipfooter=12,
        index_col=0,
        sheet_name="T17e",
        header=[0, 1, 2, 3, 4],
    ).xs("TJ", level=-1, axis=1)
    # Footnote labels lead to some strings randomly ending in numbers; we remove them here
    remove_digits = str.maketrans("", "", digits)
    # carrier names span across two column levels, which we merge with fillna
    carrier_name_func = (
        lambda x: _df.columns.to_frame()
        .iloc[:, x]
        .str.translate(remove_digits)
        .map(carriers)
    )
    _df.columns = carrier_name_func(0).fillna(carrier_name_func(1)).values

    _df = (
        _df.groupby(axis=1, level=0)
        .sum()
        .apply(utils.to_numeric)
        .rename_axis(index="year", columns="carrier_code")
    )
    _df = (
        _df.T.assign(cat_code=_df.columns.map(categories))
        .set_index("cat_code", append=True)
        .stack()
    )

    return _df


def get_ch_industry_energy_balance(path_to_excel):
    ch_subsectors = {
        1: "FC_IND_FBT_E",  # 'Food, beverages & tobacco',
        2: "FC_IND_TL_E",  # 'Textile & leather',
        3: "FC_IND_PPP_E",  # 'Paper, pulp & printing',
        4: "FC_IND_CPC_E",  # 'Chemical & petrochemical',
        5: "FC_IND_NMM_E",  # 'Non-metallic minerals',
        6: "FC_IND_NMM_E",  # 'Non-metallic minerals',
        7: "FC_IND_IS_E",  # 'Iron & steel',
        8: "FC_IND_NFM_E",  # 'Non-ferrous metals',
        9: "FC_IND_MAC_E",  # 'Machinery',
        10: "FC_IND_MAC_E",  # 'Machinery',
        11: "FC_IND_NSP_E",  # 'Not elsewhere specified (industry)', 'Wood & wood products, Mining & quarrying, Transport equipment
        12: "FC_IND_CON_E",  # 'Construction'
    }

    ch_carriers = {
        "ELEKTRIZITÄT": "E7000",  # 'electricity',
        "HEIZÖL EXTRA-LEICHT": "O4000XBIO",  # 'oil',
        "ERDGAS": "G3000",  # 'gas',
        "KOHLE": "C0000X0350-0370",  # 'solid_fuel',
        "INDUSTRIEABFÄLLE": "W6100_6220",  # 'waste',
        "HEIZÖL MITTEL UND SCHWER": "O4000XBIO",  # 'oil',
        #    'FERNWÄRME KUMULIERT': 'heat',  # total
        "FERNWÄRME BEZUG": "H8000",  # 'heat',  # purchased
        #    'FERNWÄRME ABGABE': 'heat',  # delivered
        "HOLZ": "R5110-5150_W6000RI",  # 'biofuel',
        "TOTAL": "TOTAL",
    }
    df = pd.read_excel(
        path_to_excel, sheet_name="Überblick_tot", skiprows=5, skipfooter=16, header=0
    ).dropna(how="all")

    df.index = df.index.map(df["Unnamed: 0"].where(df.TOTAL.isnull()).ffill()).map(
        ch_carriers
    )
    df = (
        df.apply(utils.to_numeric)
        .set_index("Unnamed: 0", append=True)
        .rename_axis(index=["carrier_code", "year"], columns="cat_code")
        .dropna(subset=["TOTAL"])
    )
    df.columns = (
        df.columns.str.extract("(\d+)", expand=False)
        .fillna(0)
        .astype(int)
        .map(ch_subsectors)
    )
    # combine any data that now has the same cat_code or carrier_code by using groupby
    df = df.groupby(axis=1, level=0).sum().groupby(level=[0, 1]).sum()

    return df.stack()


if __name__ == "__main__":
    generate_annual_energy_balance_nc(
        path_to_input=snakemake.input.eurostat_energy_balance,
        path_to_ch_excel=snakemake.input.ch_energy_balance,
        path_to_ch_industry_excel=snakemake.input.ch_industry_energy_balance,
        path_to_cat_names=snakemake.input.cat_names,
        path_to_carrier_names=snakemake.input.carrier_names,
        countries=snakemake.params.countries,
        path_to_result=snakemake.output[0],
    )
