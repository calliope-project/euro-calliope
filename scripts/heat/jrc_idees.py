from pathlib import Path

import numpy as np
import pandas as pd

idx = pd.IndexSlice

END_USES = {
    "Space heating": "space_heat",
    "Space cooling": "end_use_electricity",
    "Hot water": "water_heat",
    "Catering": "cooking",
}
CARRIER_NAMES = {
    "Advanced electric heating": "electricity",
    "Biomass and wastes": "biofuel",
    "Conventional electric heating": "electricity",
    "Conventional gas heaters": "gas",
    "Derived heat": "heat",
    "Electric space cooling": "electricity",
    "Electricity": "electricity",
    "Electricity in circulation and other use": "electricity",
    "Gas heat pumps": "gas",
    "Gas/Diesel oil incl. biofuels (GDO)": "oil",
    "Gases incl. biogas": "gas",
    "Geothermal energy": "renewable_heat",
    "Liquified petroleum gas (LPG)": "oil",
    "Solar": "renewable_heat",
    "Solids": "solid_fossil",
}


def process_jrc_heat_tertiary_sector_data(
    paths_to_national_data: list[str], out_path: str
):
    paths_to_national_data = [Path(p) for p in paths_to_national_data]
    dfs = []
    for file in paths_to_national_data:
        df_consumption = pd.read_excel(file, sheet_name="SER_hh_fec", index_col=0)
        df_demand = pd.read_excel(file, sheet_name="SER_hh_tes", index_col=0)
        df_summary = pd.read_excel(file, sheet_name="SER_summary", index_col=0)

        df_consumption = clean_df(df_consumption, "consumption")
        df_demand = clean_df(df_demand, "demand")

        df = pd.concat([df_consumption, df_demand])

        df_elec = (
            df_summary.loc[
                "Energy consumption by end-uses (ktoe)":"Shares of energy consumption in end-uses (in %)"
            ]
            .loc["Specific electricity uses"]
            .rename_axis(index="year")
        )

        df.loc[("electricity", "end_use_electricity"), :].update(
            df.loc[("electricity", "end_use_electricity"), :].add(df_elec, axis=1)
        )

        assert np.allclose(
            df.xs("consumption", level="energy").sum(),
            df_summary.loc[
                "Energy consumption by fuel - Eurostat structure (ktoe)"
            ].astype(float),
        )

        dfs.append(df)

    pd.concat(dfs).stack().to_csv(out_path)


def clean_df(df: pd.DataFrame, energy_type: str):
    country_code = df.index.names[0].split(" - ")[0]
    df = df.assign(end_use=np.nan)
    df.loc[df.index.isin(END_USES.keys()), "end_use"] = list(END_USES.keys())
    df.end_use = df.end_use.fillna(df.end_use.ffill())

    df = (
        df.dropna()
        .set_index("end_use", append=True)
        .drop(END_USES.keys(), level=0)
        .groupby([CARRIER_NAMES, END_USES], level=[0, 1])
        .sum()
        .assign(country_code=country_code, unit="ktoe", energy=energy_type)
        .set_index(["country_code", "unit", "energy"], append=True)
        .rename_axis(
            columns="year",
            index=["carrier_name", "end_use", "country_code", "unit", "energy"],
        )
    )
    return df


if __name__ == "__main__":
    process_jrc_heat_tertiary_sector_data(
        paths_to_national_data=snakemake.input.data, out_path=snakemake.output[0]
    )
