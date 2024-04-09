import numpy as np
import pandas as pd
from eurocalliopelib import utils

SUBSECTOR_TRANSLATIONS = {
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
    11: "FC_IND_NSP_E",  # 'Not elsewhere specified (industry)' (in the CH dataset, includes Wood & wood products, Mining & quarrying, Transport equipment)
    12: "FC_IND_CON_E",  # 'Construction'
}

ENERGY_CARRIER_TRANSLATIONS = {
    "ELEKTRIZITÄT": "E7000",  # 'electricity',
    "HEIZÖL EXTRA-LEICHT": "O4000XBIO",  # 'oil',
    "ERDGAS": "G3000",  # 'gas',
    "KOHLE": "C0000X0350-0370",  # 'solid_fuel',
    "INDUSTRIEABFÄLLE": "W6100_6220",  # 'waste',
    "HEIZÖL MITTEL UND SCHWER": "O4000XBIO",  # 'oil',
    "FERNWÄRME BEZUG": "H8000",  # 'heat',  # ASSUME: use "purchased" district heat, rather than "total" (FERNWÄRME KUMULIERT) or "delivered" (FERNWÄRME ABGABE)
    "HOLZ": "R5110-5150_W6000RI",  # 'biofuel',
    "TOTAL": "TOTAL",
}


def ch_industry_energy_balance_dataset(
    path_to_ch_industry_excel: str, path_to_output: str
):
    """
    Industry subsector energy balances are in a completely different spreadsheet,
    which are processed here.
    The spreadsheet has several tables in one sheet, with all tables sharing the
    same columns (industry subsectors) and all having years as the index.
    Each table gives values for a different energy carrier.
    We process the spreadsheet by doing the following:
    1. Extracting and renaming the per-table energy carriers (they are in the same column as the years, but have no data associated with them) and assigning those names to the index
    2. Adding years to the index as a new level
    3. Deleting all index items that relate to empty data (the rows from which the energy carrier names were extracted)
    4. Renaming the subsectors based on the subsector number

    Energy carrier and industry subsector translations are based on manual inspection.
    ASSUME: some translations are not perfect matches to eurostat energy carriers.

    """
    df = pd.read_excel(
        path_to_ch_industry_excel,
        sheet_name="Überblick_tot",
        skiprows=5,
        skipfooter=16,
        header=0,
    ).dropna(how="all")

    # 1. Extract and rename energy carrier names
    df.index = df.index.map(
        df["Unnamed: 0"].where(df.TOTAL.isnull()).ffill()
    )  # carrier names are defined over two possible columns

    # 2. Add years to index as new level and define the index/column names
    df = df.set_index("Unnamed: 0", append=True).rename_axis(
        index=["carrier_code", "year"], columns="cat_code"
    )

    # Clean data
    # commas in the data can lead to string interpretations, which we deal with here
    df = df.apply(utils.to_numeric)

    # ASSUME: we take the 'new' 2013 values from this dataset, rather than the 'old' ones since they likely account for some correction in the underlying statistical analysis.
    df = df.drop(["2013alt", "2013 alt"], level="year").rename(
        {"2013neu": 2013, "2013 neu": 2013}, level="year"
    )

    # 3. Delete all index items related to empty data
    df = df.dropna(subset=["TOTAL"])

    # Grab all-industry energy demand by carrier to test against later
    df_industry_sum = df["INDU-"].to_xarray()

    # 3. strip subsectors of all but their number (we don't use the other information)
    df.columns = df.columns.str.extract(r"(\d+)", expand=False).astype(float)

    dataarray = df.loc[:, SUBSECTOR_TRANSLATIONS.keys()].stack().to_xarray()
    dataarray = utils.rename_and_groupby(
        dataarray, SUBSECTOR_TRANSLATIONS, dim_name="cat_code"
    )
    dataarray = utils.rename_and_groupby(
        dataarray, ENERGY_CARRIER_TRANSLATIONS, dim_name="carrier_code"
    )

    # Test that we have at least captured all expected subsector energy demand
    assert np.allclose(
        dataarray.sum("cat_code", min_count=1),
        utils.rename_and_groupby(
            df_industry_sum, ENERGY_CARRIER_TRANSLATIONS, dim_name="carrier_code"
        ),
        equal_nan=True,
    )
    # Data is given in TJ, so we convert to twh here
    dataarray = utils.tj_to_twh(dataarray.assign_attrs(unit="twh"))

    dataarray.expand_dims(country_code=["CHE"]).to_netcdf(path_to_output)


if __name__ == "__main__":
    ch_industry_energy_balance_dataset(
        path_to_ch_industry_excel=snakemake.input.ch_industry_excel,
        path_to_output=snakemake.output[0],
    )
