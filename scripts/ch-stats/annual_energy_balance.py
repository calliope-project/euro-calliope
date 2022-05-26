import pandas as pd

from eurocalliopelib import utils


def ch_energy_balance(path_to_ch_energy_balance_excel, path_to_output):
    """
    Process Swiss data into a tidy dataframe that matches the structure of the annual energy balances.
    """
    ch_hh_energy_use = get_ch_energy_balance_excel_sheet(
        path_to_ch_energy_balance_excel, "T17a", skipfooter=9, cat_code="FC_OTH_HH_E"
    )
    ch_ind_energy_use = get_ch_energy_balance_excel_sheet(
        path_to_ch_energy_balance_excel, "T17b", skipfooter=12, cat_code="FC_IND_E"
    )
    ch_ser_energy_use = get_ch_energy_balance_excel_sheet(
        path_to_ch_energy_balance_excel, "T17c", skipfooter=12, cat_code="FC_OTH_CP_E"
    )

    ch_energy_use_da = utils.merge_da([
        ch_hh_energy_use,
        ch_ind_energy_use,
        ch_ser_energy_use
    ])

    ch_energy_use_da.expand_dims(country_code=["CHE"]).to_netcdf(path_to_output)


def get_ch_energy_balance_excel_sheet(path_to_excel, sheet, skipfooter, cat_code):
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


if __name__ == "__main__":
    ch_energy_balance(
        path_to_ch_energy_balance_excel=snakemake.input.ch_energy_balance_excel,
        path_to_output=snakemake.output[0]
    )
