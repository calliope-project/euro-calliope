import eurocalliopelib.utils as ec_utils
import numpy as np
import pandas as pd
import utils.jrc_idees_parser as jrc

LEVEL_ORDER = ["subsector", "country_code", "unit", "carrier", "year"]


def fill_missing_data(
    energy_balances_df: pd.DataFrame,
    cat_names_df: pd.DataFrame,
    carrier_names_df: pd.DataFrame,
    subsector_energy_consumption_df: pd.DataFrame,
    year_range: list,
):
    """
    For countries without relevant data in JRC_IDEES we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    JRC data is available until 2015, so for future years data for all countries is filled
    based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    country_codes = subsector_energy_consumption_df.index.get_level_values(
        "country_code"
    ).unique()
    if len(country_codes[0]) == 2:
        subsector_energy_consumption_df.reset_index(level="country_code", inplace=True)
        subsector_energy_consumption_df["country_code"] = (
            subsector_energy_consumption_df[
                "country_code"
            ].map(ec_utils.convert_country_code)
        )
        subsector_energy_consumption_df.set_index(
            "country_code", append=True, inplace=True
        )

    # Get annual energy balances
    industry_energy_balances = (
        energy_balances_df.unstack(["year", "country"])
        .groupby(
            [
                cat_names_df.jrc_idees.dropna().to_dict(),
                carrier_names_df.ind_carrier_name.dropna().to_dict(),
            ],
            level=["cat_code", "carrier_code"],
        )
        .sum(min_count=1)
        .sum(level="cat_code", min_count=1)
        .stack("country")
        .rename_axis(index=["cat_name", "country_code"])
        .apply(ec_utils.tj_to_ktoe)
    )

    country_codes = subsector_energy_consumption_df.index.get_level_values(
        "country_code"
    ).unique()
    # If JRC-IDEES data exists, there will be data in 'energy_consumption' for that country
    balances_with_jrc_data = industry_energy_balances[
        industry_energy_balances.index.get_level_values("country_code").isin(
            country_codes
        )
    ]
    # If JRC-IDEES data does not exist, there will only be data for that country in annual energy balances
    balances_without_jrc_data = industry_energy_balances[
        ~industry_energy_balances.index.get_level_values("country_code").isin(
            country_codes
        )
    ]

    # Compared to countries with JRC data, those without JRC-IDEES data consume X amount of energy
    # E.g. CH has ~1% of paper and pulp energy consumption compared to all countries with JRC data
    countries_without_jrc_data_contribution = balances_without_jrc_data.div(
        balances_with_jrc_data.sum(level="cat_name")
    )
    # Multiply the JRC data with the contribution from each country
    # So all JRC countries consume 3.43e4ktoe electricity in paper and pulp in 2014,
    # hence CH consumes ~3.43e2ktoe electricity in paper and pulp in 2014
    countries_without_jrc_data_consumption = (
        subsector_energy_consumption_df.sum(level=["cat_name", "unit", "carrier"])
        .align(countries_without_jrc_data_contribution)[0]
        .mul(
            subsector_energy_consumption_df.sum(
                level=["cat_name", "unit", "carrier"]
            ).align(countries_without_jrc_data_contribution)[1]
        )
    )
    all_euro_calliope_consumption = subsector_energy_consumption_df.append(
        countries_without_jrc_data_consumption.reorder_levels(
            subsector_energy_consumption_df.index.names
        )
    )
    average_consumption_per_energy_use = (
        all_euro_calliope_consumption.div(industry_energy_balances)
        .mean(axis=1)
        .where(lambda x: ~np.isinf(x) & (x > 0))
    )
    # In some instances (e.g. RO: non ferrous metals), JRC IDEES has consumption, but
    # latest Eurostat doesn't, leading to inf values
    average_consumption_per_energy_use = (
        average_consumption_per_energy_use.unstack(["cat_name", "carrier", "unit"])
        .fillna(
            average_consumption_per_energy_use.mean(
                level=["cat_name", "carrier", "unit"]
            )
        )
        .unstack()
        .reorder_levels(all_euro_calliope_consumption.index.names)
    )
    # Fill data where JRC says there is no consumption of any form in a country's industry sybsector
    # But where the energy balances show consumption (e.g. UK, Wood and wood products)
    _to_fill = all_euro_calliope_consumption.stack().unstack(["carrier", "unit"])
    _filler = (
        industry_energy_balances.stack()
        .mul(average_consumption_per_energy_use, axis=0)
        .unstack(["carrier", "unit"])
    )
    filled_jrc_no_data = _to_fill.where(_to_fill.sum(axis=1) > 0).fillna(_filler)

    # Prepare final dataframe format
    _to_fill = filled_jrc_no_data.stack(["unit", "carrier"]).unstack("year")

    # Fill data for all years where there is no JRC data
    years_to_fill = [y for y in range(*year_range) if y > jrc.MAX_YEAR]
    if years_to_fill:
        new_years_filler = _filler.stack(["unit", "carrier"]).unstack("year")
        _to_fill = _to_fill.assign(**{
            str(y): new_years_filler.loc[_to_fill.sum(axis=1) > 0, y]
            for y in years_to_fill
        })
        _to_fill.columns = _to_fill.columns.astype(int)

    # Backfill missing year data
    # (older years) / forwardfill (newer years) / linear interpolation (middle years)
    all_filled = _to_fill.interpolate(axis=1, limit_direction="both")

    return all_filled
