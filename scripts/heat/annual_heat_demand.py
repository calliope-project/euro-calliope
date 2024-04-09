from contextlib import suppress

import numpy as np
import pandas as pd
import pycountry
from eurocalliopelib import utils

END_USE_CAT_NAMES = {
    "FC_OTH_HH_E_CK": "cooking",
    "FC_OTH_HH_E_SH": "space_heat",
    "FC_OTH_HH_E_WH": "water_heat",
}

CH_ENERGY_CARRIER_TRANSLATION = {
    "Heizöl": "oil",
    "Erdgas": "gas",
    "El. Widerstandsheizungen": "direct_electric",
    "El. Wärmepumpen 1)": "heat_pump",
    "El. Ohm'sche Anlagen": "direct_electric",
    "El. Wärmepumpen": "heat_pump",
    "Elektrizität": "electricity",
    "Holz": "biofuel",
    "Kohle": "solid_fossil",
    "Fernwärme": "heat",
    "Umweltwärme": "ambient_heat",
    "Solar": "solar_thermal",
}

CH_HH_END_USE_TRANSLATION = {
    "Raumwärme": "space_heat",
    "Warmwasser": "water_heat",
    "Prozesswärme": "process_heat",
    "Beleuchtung": "end_use_electricity",
    "Klima, Lüftung, HT": "end_use_electricity",
    "I&K, Unterhaltung": "end_use_electricity",
    "Antriebe, Prozesse": "end_use_electricity",
    "sonstige": "end_use_electricity",
}

idx = pd.IndexSlice


def get_heat_demand(
    path_to_hh_end_use: str,
    path_to_ch_end_use: str,
    path_to_energy_balance: str,
    path_to_commercial_demand: str,
    path_to_carrier_names: str,
    heat_tech_params: dict[str, dict[str, float]],
    fill_missing_values: dict[str, list[str]],
    country_codes: list[str],
    path_to_electricity_demand: str,
    path_to_output: str,
) -> None:
    # Get annual energy balance data for household and commercial sectors
    energy_balance_dfs = get_energy_balances(
        path_to_energy_balance, path_to_carrier_names
    )

    # Get household final energy demand by end use
    annual_final_demand = get_household_final_energy_demand(
        path_to_hh_end_use, path_to_ch_end_use, path_to_carrier_names
    )

    # get commercial final energy demand by end use
    annual_final_demand = get_commercial_final_energy_demand(
        energy_balance_dfs["com"].add(energy_balance_dfs["oth"], fill_value=0),
        path_to_ch_end_use,
        path_to_commercial_demand,
        annual_final_demand,
        fill_missing_values=fill_missing_values,
    )

    # Fix data gaps for some countries
    annual_final_demand = hardcoded_country_cleanup(
        annual_final_demand, energy_balance_dfs["hh"]
    )

    # get electricity demand data specifically, to remove from ENTSOE timeseries
    electricity_demand = get_annual_electricity_demand(
        annual_final_demand, energy_balance_dfs
    )
    # Convert final to useful energy demand
    national_useful_heat_demand = get_national_useful_heat_demand(
        annual_final_demand, energy_balance_dfs, heat_tech_params
    )

    # Fill remaining values before saving and filter to country scope
    for df, path in zip(
        [electricity_demand, national_useful_heat_demand],
        [path_to_electricity_demand, path_to_output],
    ):
        df.stack([0, 1]).squeeze().loc[country_codes].rename("value").pipe(
            fill_remaining_missing_values
        ).to_csv(path)


def get_energy_balances(
    path_to_energy_balance: str, path_to_carrier_names: str
) -> dict[str, pd.DataFrame]:
    energy_balance_df = pd.read_csv(
        path_to_energy_balance,
        index_col=["cat_code", "carrier_code", "unit", "country", "year"],
        header=0,
    ).squeeze()
    carrier_names_df = pd.read_csv(path_to_carrier_names, index_col=0, header=0)

    balance_codes = {
        "hh": ["FC_OTH_HH_E"],
        "com": ["FC_OTH_CP_E"],
        "oth": ["FC_OTH_AF_E", "FC_OTH_FISH_E", "FC_OTH_NSP_E"],
    }
    balances = {
        sector_code: slice_energy_balance_by_sector(
            sector_code=sector_code,
            df=energy_balance_df,
            carrier_names_df=carrier_names_df,
            cat_codes=cat_codes,
        )
        for sector_code, cat_codes in balance_codes.items()
    }
    return balances


def slice_energy_balance_by_sector(
    df: pd.DataFrame,
    carrier_names_df: pd.DataFrame,
    cat_codes: list[str],
    sector_code: str,
) -> pd.DataFrame:
    df = df.loc[cat_codes]  # cat_code is always the first element

    assert df.index.get_level_values("unit").unique().tolist() == [
        "TJ"
    ], "There are other units than TJ in the energy balance data. This is not expected."

    df = (
        df.xs("TJ", level="unit")
        .apply(utils.tj_to_twh)  # TJ -> TWh
        .unstack("year")
        .rename(
            index=carrier_names_df[f"{sector_code}_carrier_name"].dropna().to_dict(),
            level="carrier_code",
        )
        .groupby(by=["carrier_code", "country"], level=["carrier_code", "country"])
        .sum()
        .loc[
            list(set(carrier_names_df[f"{sector_code}_carrier_name"].dropna()))
        ]  # remove carriers that are not relevant
        .rename_axis(["carrier_name", "country_code"], axis=0)
    )
    return df


def get_household_final_energy_demand(
    path_to_hh_end_use: str, path_to_ch_end_use: str, path_to_carrier_names: str
) -> pd.DataFrame:
    """Read data on household final energy demand."""

    carrier_names_df = pd.read_csv(path_to_carrier_names, index_col=0, header=0)

    # Index name in TSV file is 'nrg_bal,siec,unit,geo\time'
    hh_end_use_df = pd.read_csv(
        path_to_hh_end_use, delimiter="\t", index_col=0, na_values=[":", ": ", ": z"]
    )
    hh_end_use_df.index = hh_end_use_df.index.str.split(",", expand=True).rename([
        "cat_code",
        "carrier_code",
        "unit",
        "country_code",
    ])
    hh_end_use_df.columns = hh_end_use_df.columns.astype(int).rename("year")

    # remove 'countries' which are not relevant
    not_countries = [
        c
        for c in hh_end_use_df.index.get_level_values("country_code").unique()
        if len(c) > 2
    ] + ["XK"]
    hh_end_use_df = hh_end_use_df.drop(
        axis=0, level="country_code", labels=not_countries
    )

    assert check_units_removed(
        hh_end_use_df, carrier_names_df
    ), "Check that you can slice by 'TJ' only, some other units in the hh_end_use data might be relevant."

    # Just keep relevant data
    hh_end_use_df = (
        hh_end_use_df.xs("TJ", level="unit")
        .apply(utils.tj_to_twh)  # TJ -> TWh
        .dropna(how="all")
    )

    # clean up renewables info
    update_final_renewable_energy_demand(hh_end_use_df)

    country_codes_ = {
        c: utils.convert_country_code(c)
        for c in hh_end_use_df.index.get_level_values("country_code")
    }

    # Add missing renewables data to
    # rename index labels to be more readable

    hh_end_use_df = hh_end_use_df.groupby(
        [
            END_USE_CAT_NAMES,
            carrier_names_df["hh_carrier_name"].dropna().to_dict(),
            country_codes_,
        ],
        level=["cat_code", "carrier_code", "country_code"],
    ).sum()
    hh_end_use_df.index = hh_end_use_df.index.rename(
        ["end_use", "carrier_name"], level=["cat_code", "carrier_code"]
    )

    # Add Swiss data
    ch_hh_end_use_df = read_ch_hh_final_demand(path_to_ch_end_use)
    hh_end_use_df = hh_end_use_df.append(ch_hh_end_use_df, sort=True)

    # Clean up data
    hh_end_use_df = (
        hh_end_use_df.sort_index()
        .where(hh_end_use_df > 0)
        .dropna(how="all")
        .assign(cat_name="household")
        .set_index("cat_name", append=True)
    )

    return hh_end_use_df


def check_units_removed(df: pd.DataFrame, carrier_names_df: pd.DataFrame) -> bool:
    df = (  # first re-organise df
        df.stack("year")
        .unstack("unit")
        .loc[
            idx[
                END_USE_CAT_NAMES.keys(),
                carrier_names_df["hh_carrier_name"].dropna().index,
                :,
                :,
            ],
            :,
        ]
    )

    # check that when 'TJ' is NaN, the other values are also NaN.
    # Otherwise, we are missing some data. Print the df below to check.

    df = df[df["TJ"].isna() & ~(df.drop("TJ", axis=1).fillna(0) == 0).all(axis=1)]
    return len(df) == 0


def update_final_renewable_energy_demand(df: pd.DataFrame) -> None:
    """
    Some household final energy data has a higher overall renewables energy demand (RA000)
    than the sum of renewable energy carriers. Here we scale all renewable energy carriers
    evenly, to match the total of RA000.
    """

    def _get_rows_to_update(df):
        renewables = (
            df.stack()
            .unstack("carrier_code")
            .filter(regex="^R")
            .where(lambda x: x > 0)
            .dropna(how="all")
        )
        renewables_carriers = renewables.drop("RA000", axis=1).sum(axis=1, min_count=1)
        renewables_all = renewables.xs("RA000", axis=1)
        # Only update those rows where the sum of renewable energy carriers != RA000
        return renewables.loc[
            ~np.isclose(renewables_carriers, renewables_all)
        ], renewables_carriers

    to_update, renewables_carriers = _get_rows_to_update(df)
    # Some rows have no data other than RA000, so we need to assign that data to one of the
    # renewable energy carriers. We choose biofuels here (R5110-5150_W6000RI)
    completely_missing = renewables_carriers[
        renewables_carriers.isna()
    ].index.intersection(to_update.index)
    to_update.loc[completely_missing, "R5110-5150_W6000RI"] = to_update.loc[
        completely_missing, "R5110-5150_W6000RI"
    ].fillna(to_update.loc[completely_missing, "RA000"])

    # Now we scale all renewable energy carriers to match RA000
    mismatch = to_update.xs("RA000", axis=1).div(
        to_update.drop("RA000", axis=1).sum(axis=1)
    )
    updated = to_update.drop("RA000", axis=1).mul(mismatch, axis=0)
    assert np.allclose(updated.sum(axis=1), to_update.xs("RA000", axis=1))

    updated_reordered = updated.stack().unstack("year").reorder_levels(df.index.names)
    # Add new rows
    df = df.append(updated_reordered.loc[updated_reordered.index.difference(df.index)])
    # Update existing rows
    df.update(updated_reordered)
    # Ensure everything has been updated as expected
    assert _get_rows_to_update(df)[0].empty


def read_ch_hh_final_demand(path_to_ch_end_use: str) -> pd.DataFrame:
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """
    space_heat = get_ch_sheet(
        path_to_ch_end_use,
        "Tabelle 18",
        skipfooter=8,
        translation=CH_ENERGY_CARRIER_TRANSLATION,
    )
    water_heat = get_ch_sheet(
        path_to_ch_end_use,
        "Tabelle 20",
        skipfooter=5,
        translation=CH_ENERGY_CARRIER_TRANSLATION,
    )
    # Quirk of the excel is that there is no space in this sheet name
    cooking = get_ch_sheet(
        path_to_ch_end_use,
        "Tabelle21",
        skipfooter=4,
        translation=CH_ENERGY_CARRIER_TRANSLATION,
    )

    df = (
        pd.concat(
            [space_heat, water_heat, cooking],
            keys=("space_heat", "water_heat", "cooking"),
            names=["cat_name", "carrier_name"],
        )
        .assign(country_code="CHE")
        .set_index("country_code", append=True)
    )

    # Columns are years
    df.columns = df.columns.astype(int).rename("year")

    return df.apply(utils.pj_to_twh)  # PJ -> TWh


def get_commercial_final_energy_demand(
    energy_balance: pd.DataFrame,
    path_to_ch_end_use: str,
    path_to_jrc_end_use: str,
    annual_final_energy_demand: pd.DataFrame,
    fill_missing_values: dict[str, list[str]],
) -> pd.DataFrame:
    """
    Use JRC IDEES service sector final energy demand to estimate demand for
    space heating and hot water in the commercial sector across all countries.
    Add Swiss data from Swiss govt. stats.
    """

    # 'fuel' is actually just generic non-electric energy, which distribute based on household data
    ch_con_fuel = read_ch_non_hh_non_electricity_demand(
        path_to_ch_end_use, "Tabelle 25", annual_final_energy_demand
    )
    ch_con_elec = read_ch_non_hh_electricity_demand(path_to_ch_end_use, "Tabelle26")

    jrc_end_use_df = (
        pd.read_csv(
            path_to_jrc_end_use,
            index_col=[
                "carrier_name",
                "end_use",
                "country_code",
                "unit",
                "energy",
                "year",
            ],
        )
        .squeeze()
        .xs("final_energy", level="energy")
    )

    # Map JRC end uses to annual commercial demand
    mapped_end_uses = map_jrc_to_eurostat(
        energy_balance, jrc_end_use_df, fill_missing_values
    )

    # Add Swiss data and ambient heat from heat pumps
    mapped_end_uses = (
        mapped_end_uses.append(
            ch_con_fuel.rename({"process_heat": "cooking"}).reorder_levels(
                mapped_end_uses.index.names
            )
        )
        .append(
            ch_con_elec.rename({"process_heat": "cooking"}).reorder_levels(
                mapped_end_uses.index.names
            )
        )
        .append(
            energy_balance.loc[  # JRC data only refers to heat pumps for heating in space heating
                ["ambient_heat"]
            ]
            .assign(end_use="space_heat")
            .set_index("end_use", append=True)
            .stack()
            .reorder_levels(mapped_end_uses.index.names)
        )
    )
    mapped_end_uses.index = mapped_end_uses.index.remove_unused_levels()

    annual_final_energy_demand = annual_final_energy_demand.append(
        mapped_end_uses.where(mapped_end_uses > 0)
        .dropna()
        .unstack("year")[annual_final_energy_demand.columns]
        .assign(cat_name="commercial")
        .set_index("cat_name", append=True)
        .reorder_levels(annual_final_energy_demand.index.names)
    ).sort_index()

    return annual_final_energy_demand


def map_jrc_to_eurostat(
    energy_balance: pd.DataFrame,
    jrc_end_use_df: pd.DataFrame,
    fill_missing_values: dict[str, list[str]],
) -> pd.DataFrame:
    jrc_end_use_df = (
        jrc_end_use_df.xs("ktoe", level="unit")
        .rename(utils.convert_country_code, level="country_code")
        .apply(utils.ktoe_to_twh)  # kTOE -> TWh
    )
    jrc_end_use_percent = (
        jrc_end_use_df.div(jrc_end_use_df.unstack("end_use").sum(axis=1))
        .unstack("year")
        .dropna(how="all")
        .stack("year")
    )

    jrc_end_use_percent = fill_missing_countries_and_years(
        jrc_data=jrc_end_use_percent, fill_missing_values=fill_missing_values
    )

    mapped_end_uses = (
        jrc_end_use_percent.align(energy_balance.stack())[1]
        .mul(jrc_end_use_percent)
        .dropna()
    )

    return mapped_end_uses


def fill_missing_countries_and_years(
    jrc_data: pd.DataFrame, fill_missing_values: dict[str, str]
) -> pd.DataFrame:
    fill_missing_values.pop(
        "CHE", None
    )  # do not fill Swiss values as they are taken from CH data
    jrc_data = jrc_data.unstack("country_code")
    with suppress(
        KeyError
    ):  # it's fine. Just checking there is no MultiIndex in the columns
        jrc_data = jrc_data.loc[:, "value"]
    for country, neighbors in fill_missing_values.items():
        jrc_data = jrc_data.assign(**{country: jrc_data[neighbors].mean(axis=1)})

    jrc_data = jrc_data.stack().unstack("year")
    jrc_data = jrc_data.assign(**{
        str(year): jrc_data[2015] for year in range(2016, 2019)
    })
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()


def read_ch_non_hh_electricity_demand(
    path_to_ch_end_use: str,
    sheet_name: str,
) -> pd.DataFrame:
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """

    return (
        get_ch_sheet(
            path_to_ch_end_use,
            sheet_name,
            skipfooter=4,
            translation=CH_HH_END_USE_TRANSLATION,
        )
        .assign(carrier_name="electricity", country_code="CHE")
        .set_index(["country_code", "carrier_name"], append=True)
        .stack()
        .rename_axis(index=["end_use", "country_code", "carrier_name", "year"])
        .apply(utils.pj_to_twh)  # PJ -> TWh
    )


def read_ch_non_hh_non_electricity_demand(
    path_to_ch_end_use: str, sheet_name: str, hh_final_energy_demand: pd.DataFrame
):
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """
    ch_con = get_ch_sheet(
        path_to_ch_end_use,
        sheet_name,
        skipfooter=4,
        translation=CH_HH_END_USE_TRANSLATION,
    )
    # this is actually just generic non-electric energy,
    # which we assign to fuels using household ratios
    # ASSUME Swiss carrier ratios in commerce are the same as in households
    hh_con = hh_final_energy_demand.xs(
        ("CHE", "household"), level=("country_code", "cat_name")
    )
    hh_ratios = hh_con.div(
        hh_con.drop("electricity", level="carrier_name").sum(level=["end_use"]),
        axis=0,
    )
    ch_con_disaggregated = hh_ratios.mul(ch_con, level="end_use", axis=0).dropna(
        how="all"
    )

    assert np.allclose(ch_con_disaggregated.sum(level="end_use"), ch_con)
    ch_con = ch_con_disaggregated.reset_index("carrier_name")
    return (
        ch_con.assign(country_code="CHE")
        .set_index(["country_code", "carrier_name"], append=True)
        .stack()
        .rename_axis(index=["end_use", "country_code", "carrier_name", "year"])
        .apply(utils.pj_to_twh)  # PJ -> TWh
    )


def hardcoded_country_cleanup(
    annual_final_energy_demand: pd.DataFrame, energy_balance_df: pd.DataFrame
) -> pd.DataFrame:
    # Montenegro data isn't in the database, so we combine estimates from the World bank
    # with average end use contributions of neighbouring countries
    # biomass fuel = 69%, electricity = 28%, oil/solid fossil fuel = 1-2%
    # See The World Bank Montenegro Second Energy Efficiency Project (P165509)
    # FIXME this function is surprising. Why is a special case treatment necessary for Montenegro?
    # Montenegro is missing from JRC IDEES like all Western Balkan countries.
    # The data that is supposed to be filled in this function is in fact existing.
    # See https://github.com/calliope-project/euro-calliope/issues/298.
    MNE_energy_balance = energy_balance_df.loc[idx[:, "MNE"], :]
    MNE_heat_electricity_demand = MNE_energy_balance.loc["biofuel"] * 0.28 / 0.69
    MNE_energy_balance.loc["electricity"].update(
        MNE_heat_electricity_demand
    )  # FIXME this line is broken. It does not do anything.
    neighbours = ["SRB", "HRV", "ALB", "BIH"]
    neighbour_demand = annual_final_energy_demand.loc[
        idx[END_USE_CAT_NAMES.values(), :, neighbours, ["household"]], :
    ]
    end_use_contributions = (
        neighbour_demand.sum(level=["end_use", "country_code", "cat_name"])
        .div(neighbour_demand.sum(level="country_code"))
        .mean(level=["end_use", "cat_name"])
    )
    MNE_end_use = (
        end_use_contributions.stack()
        .unstack(["end_use", "cat_name"])
        .mul(MNE_energy_balance.stack(), axis=0)
        .stack([0, 1])
        .unstack("year")
        .reorder_levels(annual_final_energy_demand.index.names)
    )
    MNE_end_use = MNE_end_use.where(MNE_end_use > 0).dropna(how="all")
    annual_final_energy_demand = annual_final_energy_demand.append(
        MNE_end_use, sort=True
    ).sort_index()

    return annual_final_energy_demand


def fill_data_gaps(
    end_use_df: pd.DataFrame, energy_balance_dfs: dict[str, pd.DataFrame], fill: str
) -> pd.DataFrame:
    end_use_df = end_use_df.where(end_use_df > 0)

    # Fill the household sector's end-use data based on total sectoral demand
    hh_country_energy_balance = (
        energy_balance_dfs["hh"]
        .sum(level="country_code")
        .stack()
        .where(lambda x: x > 0)
    )

    end_use_df.loc[:, idx[:, "household"]] = end_use_df.loc[
        :, idx[:, "household"]
    ].fillna(
        end_use_df.loc[:, idx[:, "household"]]
        .div(hh_country_energy_balance, axis=0)
        .mean(level="country_code")
        .mul(hh_country_energy_balance, level="country_code", axis=0)
    )

    end_use_df = end_use_df.where(end_use_df > 0)

    # For all remaining gaps, fill with mean/first/last/max/min/whatever data for the
    # country, based on the string given by 'fill'
    end_use_df = end_use_df.fillna(end_use_df.groupby("country_code").agg(fill))

    return end_use_df


def get_annual_electricity_demand(
    annual_final_energy_demand: pd.DataFrame,
    energy_balance_dfs: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    """
    Get annual energy demand coming from electricity, to remove later from the
    electricity demand profile.
    Gaps in electricity demand are filled before saving, based on annual energy
    demand
    """
    electricity_demand = (
        annual_final_energy_demand.drop("end_use_electricity")
        .loc[idx[:, ["electricity", "direct_electric", "heat_pump"], :, :], :]
        .sum(level=["end_use", "country_code", "cat_name"])
        .stack()
        .unstack(["end_use", "cat_name"])
    )

    electricity_demand = fill_data_gaps(
        electricity_demand, energy_balance_dfs, fill="first"
    )
    electricity_demand = electricity_demand.rename(
        lambda x: f"{x}_historically_electrified", level="end_use", axis=1
    )
    return electricity_demand


def get_national_useful_heat_demand(
    annual_final_energy_demand: pd.DataFrame,
    energy_balance_dfs: dict[str, pd.DataFrame],
    heat_tech_params: dict[str, dict[str, float]],
) -> pd.DataFrame:
    """
    Derive useful heat demand from final energy demand.
    """

    demands = []
    for end_use in ["space_heat", "water_heat", "cooking"]:
        _demand = (
            annual_final_energy_demand.loc[[end_use]]
            .mul(efficiencies(heat_tech_params[end_use]), level="carrier_name", axis=0)
            .sum(level=["end_use", "country_code", "cat_name"])
        )
        # sense check: useful demand is at least the minimum efficiency * final demand
        assert (
            (
                _demand
                >= (
                    annual_final_energy_demand.loc[[end_use]]
                    .mul(efficiencies(heat_tech_params[end_use]).min())
                    .sum(level=["end_use", "country_code", "cat_name"])
                )
            )
            .all()
            .all()
        )
        demands.append(_demand)

    demand = pd.concat(demands).stack().unstack(["end_use", "cat_name"])

    # Fill gaps in demand using energy balance data
    # This is done by finding the average contribution of cooking, space heating
    # and water heating to HH demand (from energy balances),
    # then applying to years in which HH sub class data is not available.
    demand = fill_data_gaps(demand, energy_balance_dfs, fill="mean")

    return demand


def efficiencies(params: dict[str, float]) -> pd.Series:
    return pd.Series({
        "biogas": params.get("gas-eff", np.nan),
        "biofuel": params.get("biofuel-eff", np.nan),
        "solid_fossil": params.get("solid-fossil-eff", np.nan),
        "natural_gas": params.get("gas-eff", np.nan),
        "manufactured_gas": params.get("gas-eff", np.nan),
        "gas": params.get("gas-eff", np.nan),
        "oil": params.get("oil-eff", np.nan),
        "solar_thermal": params.get("solar-thermal-eff", np.nan),
        "renewable_heat": params.get("solar-thermal-eff", np.nan),
        "electricity": params.get("electricity-eff", np.nan),
        "direct_electric": 1,  # don't need to deal with heat pump COP if direct electric is 100% efficient
        "heat": 1,
        # heat demand met by heat pumps = heat pump electricity + ambient heat
        "heat_pump": 1,
        "ambient_heat": 1,
    })


def get_ch_sheet(
    path_to_excel: str, sheet: str, skipfooter, translation=None
) -> pd.DataFrame:
    df = pd.read_excel(
        path_to_excel, sheet_name=sheet, skiprows=9, skipfooter=skipfooter, index_col=1
    ).drop(["Unnamed: 0", "Δ ’00 – ’18"], axis=1)
    df.index = df.index.str.strip()
    df.columns = df.columns.astype(int)
    df = df.drop(2019, axis=1, errors="ignore")

    if translation is not None:
        return df.groupby(translation).sum()
    else:
        return df


def fill_remaining_missing_values(
    df: pd.DataFrame,
) -> pd.DataFrame:
    # FIXME we shouldn't fill NaNs with 0 because this means 0 demand, maybe this should be in pre-processing
    return df.fillna(0)


if __name__ == "__main__":
    get_heat_demand(
        path_to_hh_end_use=snakemake.input.hh_end_use,
        path_to_ch_end_use=snakemake.input.ch_end_use,
        path_to_energy_balance=snakemake.input.energy_balance,
        path_to_commercial_demand=snakemake.input.commercial_demand,
        path_to_carrier_names=snakemake.input.carrier_names,
        heat_tech_params=snakemake.params.heat_tech_params,
        country_codes=[
            pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
        ],
        fill_missing_values=snakemake.params.fill_missing_values,
        path_to_electricity_demand=snakemake.output.electricity,
        path_to_output=snakemake.output.total_demand,
    )
