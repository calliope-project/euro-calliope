import pandas as pd

from eurocalliopelib import utils

idx = pd.IndexSlice
YEARS = range(2000, 2019)


def generate_annual_energy_balance_nc(
    path_to_input,
    path_to_cat_names,
    path_to_carrier_names,
    path_to_ch_excel,
    path_to_ch_industry_excel,
    countries,
    path_to_result,
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

    index_names = ["cat_code", "carrier_code", "unit", "country_code"]
    da = utils.read_eurostat_tsv(path_to_input, index_names).stack().to_xarray()

    da = da.loc[{
        "cat_code": cat_names.index.intersection(da.cat_code),
        "carrier_code": carrier_names.index.intersection(da.carrier_code),
        "unit": "TJ"
    }]

    country_code_mapping = valid_countries_to_alpha3(da.country_code)
    da = utils.rename_and_groupby(da, country_code_mapping, dim="country_code")

    da = utils.tj_to_twh(da).drop_vars("unit").assign_attrs({"unit": "twh"})

    # Add CH energy use (only covers a subset of sectors and carriers, but should be enough)
    ch_energy_use_da = add_ch_energy_balance(path_to_ch_excel, path_to_ch_industry_excel)

    all_da = (
        utils.merge_da([da, ch_energy_use_da], "annual_energy_balances")
        .sel(year=YEARS)
    )

    # Ensure we have at least the countries in the model scope defined.
    # All the others are worth keeping to do data gap filling later on in the workflow.
    assert all(utils.convert_country_code(country) in da.country_code.to_index() for country in countries)

    all_da.to_netcdf(path_to_result)


def add_ch_energy_balance(path_to_ch_excel, path_to_ch_industry_excel):
    """
    Process Swiss data into a tidy dataframe that matches the structure of the annual energy balances.
    """
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

    ch_energy_use_da = utils.merge_da([
        ch_hh_energy_use,
        ch_ind_energy_use,
        ch_ser_energy_use,
        ch_waste_energy_use,
        ch_industry_subsector_energy_use,
        ch_transport_energy_use,
    ])

    return ch_energy_use_da.expand_dims(country_code=["CHE"])


def get_ch_energy_balance_sheet(path_to_excel, sheet, skipfooter, cat_code):
    """
    Get energy balance data per sector, which requires translating energy carrier
    names and removing footnotes that are in the spreadsheet

    Parameters
    ----------
    sheet : sheet name in the excel
    skipfooter : number of footer rows to ignore which correspond to footnotes in the sheet
    cat_code : sector (category) for which the sheet has data, corresponding to eurostat sectors.
    """
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
    df = (
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
    df.columns = (
        df
        .columns
        .get_level_values(0)
        .str
        .translate(utils.remove_digits())
        .map(ch_energy_carriers)
        .rename("carrier_code")
    )
    df.index.rename("year", inplace=True)

    df_twh = (
        df
        .apply(utils.to_numeric)
        .apply(utils.tj_to_twh)
        .stack()
    )

    da = df_twh.to_xarray().expand_dims(cat_code=[cat_code]).assign_attrs({"unit": "twh"})

    return da


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
    waste_stream_twh = (
        waste_stream_gwh
        .rename_axis(index="year")
        .apply(utils.gwh_to_twh)
    )
    waste_stream_da = (
        waste_stream_twh
        .to_xarray()
        .expand_dims(carrier_code=["W6100_6220"], cat_code=["TI_EHG_E"])
    )
    return waste_stream_da.assign_attrs({"unit": "twh"})


def get_ch_transport_energy_balance(path_to_excel):
    """
    Swiss transport sector energy balance sheet is structured differently to the other sectors, requiring inference of what end uses the fuels are for (road, rail, and aviation).
    ASSUME: petrol, diesel, and gas for road transport.
    ASSUME: electricity for rail.
    ASSUME: kerosene for aviation.
    """
    carriers_categories = {
        "Elektrizität": ("E7000", "FC_TRA_RAIL_E"),
        "Gas übriger Vekehr": ("G3000", "FC_TRA_ROAD_E"),
        "Übrige erneuerbare Energien": ("R5220B", "FC_TRA_ROAD_E"),
        "davon Benzin": ("O4652XR5210B", "FC_TRA_ROAD_E"),
        "davon Diesel": ("O4671XR5220B", "FC_TRA_ROAD_E"),
        "davon Flugtreibstoffe": ("O4000XBIO", "INTAVI"),
    }
    df = pd.read_excel(
        path_to_excel,
        skiprows=6,
        skipfooter=12,
        index_col=0,
        sheet_name="T17e",
        header=[0, 1, 2, 3, 4],
    ).xs("TJ", level=-1, axis=1)
    # Footnote labels lead to some strings randomly ending in numbers; we remove them here
    # carrier names span across two column levels, which we merge with fillna
    carrier_name_func = (
        lambda x: df.columns.to_frame()
        .iloc[:, x]
        .str.translate(utils.remove_digits())
        .map(carriers_categories)
    )
    df.columns = carrier_name_func(0).fillna(carrier_name_func(1)).values

    df_twh = (
        df
        .groupby(axis=1, level=0).sum()
        .apply(utils.to_numeric)
        .apply(utils.tj_to_twh)
        .rename_axis(index="year")
    )
    df_twh.columns = pd.MultiIndex.from_tuples(df_twh.columns).rename(["carrier_code", "cat_code"])

    da = df_twh.stack(["carrier_code", "cat_code"]).to_xarray()

    return da.assign_attrs({"unit": "twh"})


def get_ch_industry_energy_balance(path_to_excel):
    """
    Industry subsector energy balances are in a completely different spreadsheet,
    which are processed here.
    ASSUME: some translations are not perfect matches to eurostat energy carriers.
    """
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
        # 'Not elsewhere specified (industry)', 'Wood & wood products, Mining & quarrying, Transport equipment
        11: "FC_IND_NSP_E",
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
    df_twh = (
        df
        .set_index("Unnamed: 0", append=True)
        .apply(utils.to_numeric)
        .apply(utils.tj_to_twh)
        .rename_axis(index=["carrier_code", "year"], columns="cat_code")
        # ASSUME: we take the 'new' 2013 values from this dataset, rather than the 'old' ones
        # since they likely account for some correction in the underlying statistical analysis.
        .drop(["2013alt", "2013 alt"], level="year")
        .rename({"2013neu": 2013, "2013 neu": 2013}, level="year")
        .dropna(subset=["TOTAL"])
    )
    df_twh.columns = (
        df_twh
        .columns
        .str
        .extract("(\d+)", expand=False)
        .fillna(0)
        .astype(int)
        .map(ch_subsectors)
    )
    # combine any data that now has the same cat_code or carrier_code by using groupby
    df_twh = (
        df_twh
        .groupby(axis=1, level="cat_code").sum()
        .groupby(level=["carrier_code", "year"]).sum()
    )

    da = df_twh.stack().to_xarray()

    return da.assign_attrs({"unit": "twh"})


def valid_countries_to_alpha3(country_codes):
    mapped_codes = {}
    for country_code in country_codes.values:
        try:
            mapped_codes[country_code] = utils.convert_country_code(country_code)
        except LookupError:
            print(f"Skipping country/region {country_code} in annual energy balances")
            continue
    return mapped_codes


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
