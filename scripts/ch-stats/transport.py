import pandas as pd
from eurocalliopelib import utils

CARRIER_TO_CARRIER_CODE_AND_END_USE = {
    "Elektrizität": ("E7000", "FC_TRA_RAIL_E"),
    "Gas übriger Vekehr": ("G3000", "FC_TRA_ROAD_E"),
    "Übrige erneuerbare Energien": ("R5220B", "FC_TRA_ROAD_E"),
    "davon Benzin": ("O4652XR5210B", "FC_TRA_ROAD_E"),
    "davon Diesel": ("O4671XR5220B", "FC_TRA_ROAD_E"),
    "davon Flugtreibstoffe": ("O4000XBIO", "INTAVI"),
}


def ch_transport_energy_balance(path_to_ch_energy_balance_excel, path_to_output):
    """
    Swiss transport sector energy balance sheet is structured differently to the other sectors,
    requiring us to infer the end uses for the fuels (road, rail, and aviation).

    ASSUME: petrol, diesel, and gas for road transport.
    ASSUME: electricity for rail.
    ASSUME: kerosene for international aviation.
    """
    df = pd.read_excel(
        path_to_ch_energy_balance_excel,
        skiprows=6,
        skipfooter=12,
        index_col=0,
        sheet_name="T17e",
        header=[0, 1, 2, 3, 4],
    ).xs("TJ", level=-1, axis=1)

    # carrier names span two column levels, which we merge with fillna
    carrier_name_func = (
        lambda x: df.columns.to_frame()
        .iloc[:, x]
        .str.translate(
            utils.remove_digits()
        )  # Footnote labels lead to some strings randomly ending in numbers; we remove them here
        .map(CARRIER_TO_CARRIER_CODE_AND_END_USE)
    )
    df.columns = carrier_name_func(0).fillna(carrier_name_func(1)).values

    df_twh = (
        df.groupby(axis=1, level=0)
        .sum()
        .apply(utils.to_numeric)
        .apply(utils.tj_to_twh)
        .rename_axis(index="year")
    )
    df_twh.columns = pd.MultiIndex.from_tuples(df_twh.columns).rename([
        "carrier_code",
        "cat_code",
    ])

    da = df_twh.stack(["carrier_code", "cat_code"]).to_xarray()

    da.expand_dims(country_code=["CHE"]).assign_attrs(unit="twh").to_netcdf(
        path_to_output
    )


if __name__ == "__main__":
    ch_transport_energy_balance(
        path_to_ch_energy_balance_excel=snakemake.input.ch_energy_balance_excel,
        path_to_output=snakemake.output[0],
    )
