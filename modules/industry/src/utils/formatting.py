import numpy as np
import pandas as pd
import eurocalliopelib.utils as utils


def fill_missing_data(energy_balances, cat_names, carrier_names, energy_consumption):
    """
    There are 7 countries without relevant data in JRC_IDEES, so we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    2016-2018 data for all 35 countries is filled in based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    country_codes = energy_consumption.index.get_level_values("country_code").unique()
    if len(country_codes[0]) == 2:
        energy_consumption.reset_index(level="country_code", inplace=True)
        energy_consumption["country_code"] = energy_consumption["country_code"].map(utils.convert_country_code)
        energy_consumption.set_index("country_code", append=True, inplace=True)

    # Get annual energy balances
    industry_energy_balances = (
        energy_balances
        .unstack(['year', 'country'])
        .groupby([
            cat_names.jrc_idees.dropna().to_dict(),
            carrier_names.ind_carrier_name.dropna().to_dict()
        ], level=['cat_code', 'carrier_code']).sum(min_count=1)
        .sum(level='cat_code', min_count=1)
        .stack('country')
        .rename_axis(index=['cat_name', 'country_code'])
        .apply(utils.tj_to_ktoe)
    )

    # If JRC-IDEES data exists, there will be data in 'energy_consumption' for that country
    jrc_country_codes = energy_consumption.index.get_level_values("country_code").unique()
    balances_with_jrc_data = industry_energy_balances[industry_energy_balances.index.get_level_values("country_code").isin(jrc_country_codes)]
    # If JRC-IDEES data does not exist, there will only be data for that country in annual energy balances
    balances_without_jrc_data = industry_energy_balances[~industry_energy_balances.index.get_level_values("country_code").isin(jrc_country_codes)]

    # Compared to countries with JRC data, those without JRC-IDEES data consume X amount of energy
    # E.g. CH has ~1% of paper and pulp energy consumption compared to all countries with JRC data
    countries_without_jrc_data_contribution = balances_without_jrc_data.div(
        balances_with_jrc_data.sum(level='cat_name')
    )
    # Multiply the JRC data with the contribution from each country
    # So all JRC countries consume 3.43e4ktoe electricity in paper and pulp in 2014,
    # hence CH consumes ~3.43e2ktoe electricity in paper and pulp in 2014
    countries_without_jrc_data_consumption = (
        energy_consumption
        .sum(level=['cat_name', 'unit', 'carrier'])
        .align(countries_without_jrc_data_contribution)[0]
        .mul(energy_consumption
             .sum(level=['cat_name', 'unit', 'carrier'])
             .align(countries_without_jrc_data_contribution)[1])
    )
    breakpoint()
    all_euro_calliope_consumption = energy_consumption.append(
        countries_without_jrc_data_consumption
        .reorder_levels(energy_consumption.index.names)
    )
    average_consumption_per_energy_use = (
        all_euro_calliope_consumption
        .div(industry_energy_balances)
        .mean(axis=1)
        .where(lambda x: ~np.isinf(x) & (x > 0))
    )
    # In some instances (e.g. RO: non ferrous metals), JRC IDEES has consumption, but
    # latest Eurostat doesn't, leading to inf values
    average_consumption_per_energy_use = (
        average_consumption_per_energy_use
        .unstack(['cat_name', 'carrier', 'unit'])
        .fillna(
            average_consumption_per_energy_use
            .mean(level=['cat_name', 'carrier', 'unit'])
        )
        .unstack()
        .reorder_levels(all_euro_calliope_consumption.index.names)
    )
    # Fill data where JRC says there is no consumption of any form in a country's industry sybsector
    # But where the energy balances show consumption (e.g. UK, Wood and wood products)
    breakpoint()
    _to_fill = all_euro_calliope_consumption.stack().unstack(["carrier", "unit"])
    _filler = industry_energy_balances.stack().mul(average_consumption_per_energy_use, axis=0).unstack(["carrier", "unit"])
    filled_jrc_no_data = _to_fill.where(_to_fill.sum(axis=1) > 0).fillna(_filler)

    # Fill data for all years between 2016 and 2018, since there is no data from JRC
    _to_fill = filled_jrc_no_data.stack(["unit", "carrier"]).unstack("year")
    new_years_filler = _filler.stack(["unit", "carrier"]).unstack("year")
    filled_2016_to_2018 = _to_fill.assign(**{str(_year): new_years_filler.loc[_to_fill.sum(axis=1) > 0, _year] for _year in range(2016, 2019)})
    filled_2016_to_2018.columns = filled_2016_to_2018.columns.astype(int)

    # Fill any remaining missing data (e.g. BA, pre 2013) with backfill (older years) / forwardfill (newer years) / linear interpolation (middle years)
    all_filled = filled_2016_to_2018.interpolate(axis=1, limit_direction="both")

    return all_filled


def verify_data(euro_calliope_consumption, industry_energy_balances, year_range):
    """
    Check that all our processing of data leads to sensible results relative to
    industry energy balances. What we mean by this is that > 95% of total euro-calliope
    energy demand should be within 0.5x and 1.5x historical primary energy
    consumption across subsectors. It's a crude estimate, but should catch cases
    where we're completely off / have missed some energy demand.
    """
    # Get the ratio of euro-calliope energy consumption to eurostat energy consumption
    consumption_diff = (
        euro_calliope_consumption
        .xs('ktoe', level='unit')
        .sum(level=['cat_name', 'country_code'], min_count=1)
        .div(industry_energy_balances.where(industry_energy_balances > 0))
        .loc[:, year_range]
        .stack()
    )
    # If industry energy balances have finite values, so does euro-calliope
    assert consumption_diff[consumption_diff == 0].sum() == 0

    # Calculate contribution of outliers in the consumption ratios
    Q1 = consumption_diff.quantile(0.25)
    Q3 = consumption_diff.quantile(0.75)
    IQR = Q3 - Q1
    # We would consider ratios outside 0.5-1.5 to be outliers, so just make sure the
    # threshold for outliers in the data is at least better than that
    assert Q1 - IQR * 3 >= 0.5
    assert Q3 + IQR * 3 <= 1.5
    outliers = consumption_diff[(consumption_diff < 0.5) | (consumption_diff > 1.5)]

    # We're happy with < 10% outliers
    assert len(outliers) / len(consumption_diff) <= 0.10

    # We're happy with < 5% of total demand being within the outliers
    assert (
        industry_energy_balances.stack().reindex(outliers.index).sum() /
        industry_energy_balances.loc[:, year_range].stack().sum() <= 0.05
    )
