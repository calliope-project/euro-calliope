import pandas as pd
import numpy as np

from eurocalliopelib import utils

END_USE_CAT_NAMES = {
    'FC_OTH_HH_E_CK': 'cooking',
    'FC_OTH_HH_E_SH': 'space_heat',
    'FC_OTH_HH_E_WH': 'water_heat'
}

CH_ENERGY_CARRIER_TRANSLATION = {
    'Heizöl': 'oil',
    'Erdgas': 'gas',
    'El. Widerstandsheizungen': 'direct_electric',
    'El. Wärmepumpen 1)': 'heat_pump',
    "El. Ohm'sche Anlagen": 'direct_electric',
    'El. Wärmepumpen': 'heat_pump',
    'Elektrizität': 'electricity',
    'Holz': 'biofuel',
    'Kohle': 'solid_fossil',
    'Fernwärme': 'heat',
    'Umweltwärme': 'ambient_heat',
    'Solar': 'solar_thermal',
}

CH_HH_END_USE_TRANSLATION = {
    'Raumwärme': 'space_heat',
    'Warmwasser': 'water_heat',
    'Prozesswärme': 'process_heat',
    'Beleuchtung': 'end_use_electricity',
    'Klima, Lüftung, HT': 'end_use_electricity',
    'I&K, Unterhaltung': 'end_use_electricity',
    'Antriebe, Prozesse': 'end_use_electricity',
    'sonstige': 'end_use_electricity'
}

idx = pd.IndexSlice


def get_heat_demand(
    path_to_hh_end_use, path_to_ch_end_use, path_to_energy_balance, path_to_commercial_demand,
    path_to_carrier_names, countries, heat_tech_params,
    path_to_electricity_consumption, path_to_output
):

    country_codes = {
        utils.get_alpha2(i, eurostat=True): utils.get_alpha3(i) for i in countries
    }
    # Get annual energy balance data for household and commercial sectors
    energy_balance_dfs = get_energy_balances(
        path_to_energy_balance, path_to_carrier_names, country_codes
    )

    # Get household energy consumption by end use
    annual_consumption = get_household_energy_consumption(
        path_to_hh_end_use, path_to_ch_end_use, path_to_carrier_names, country_codes
    )

    # get commercial energy consumption by end use
    annual_consumption = get_commercial_energy_consumption(
        energy_balance_dfs['com'].add(energy_balance_dfs['oth'], fill_value=0),
        path_to_ch_end_use, path_to_commercial_demand, country_codes, annual_consumption
    )

    # Fix data gaps for some countries
    annual_consumption = hardcoded_country_cleanup(
        annual_consumption, energy_balance_dfs['hh']
    )

    # get electricity consumption data specifically, to remove from ENTSOE timeseries
    electricity_consumption = get_annual_electricity_consumption(
        annual_consumption, energy_balance_dfs
    )
    # Convert consumption to demand for heat
    national_heat_demand = get_national_heat_demand(
        annual_consumption, energy_balance_dfs, heat_tech_params
    )

    # Set unit before saving
    def set_unit(df, unit='twh'):
        return df.assign(unit=unit).set_index('unit', append=True)

    set_unit(electricity_consumption).stack([0, 1]).to_csv(path_to_electricity_consumption)
    set_unit(national_heat_demand).stack([0, 1]).to_csv(path_to_output)


def get_energy_balances(energy_balance, carrier_names, country_codes):
    energy_balance_df = utils.read_tdf(energy_balance)
    carrier_names_df = pd.read_csv(carrier_names, index_col=0, header=0)

    balances = {}
    balance_codes = {
        'hh': ['FC_OTH_HH_E'],
        'com': ['FC_OTH_CP_E'],
        'oth': ['FC_OTH_AF_E', 'FC_OTH_FISH_E', 'FC_OTH_NSP_E']
    }
    for balance_code, cat_code in balance_codes.items():
        _df = energy_balance_df.loc[cat_code]  # cat_code is always the first element
        _df = (
            _df.xs('TJ', level='unit')
            .apply(utils.tj_to_twh)  # TJ -> TWh
            .unstack('year')
            .groupby(
                [carrier_names_df[f'{balance_code}_carrier_name'].dropna().to_dict(), country_codes],
                level=['carrier_code', 'country']
            )
            .sum()
            .rename_axis(['carrier_name', 'country_code'], axis=0)
        )
        balances[balance_code] = _df

    return balances


def get_household_energy_consumption(
    hh_end_use, ch_end_use, carrier_names, country_codes
):
    """
    Combine annual energy balance information with data on household
    end use energy consumption
    """

    carrier_names_df = pd.read_csv(carrier_names, index_col=0, header=0)

    # Index name in TSV file is 'nrg_bal,siec,unit,geo\time'
    hh_end_use_df = utils.read_eurostat_tsv(
        hh_end_use, ['cat_code', 'carrier_code', 'unit', 'country_code']
    )

    # Just keep relevant data
    hh_end_use_df = (
        hh_end_use_df.xs('TJ', level='unit')
        .apply(utils.tj_to_twh)  # TJ -> TWh
        .astype(float)
        .dropna(how='all')
    )
    # clean up renewables info
    update_renewable_energy_consumption(hh_end_use_df)

    # Add missing renewables data to
    # rename index labels to be more readable
    hh_end_use_df = hh_end_use_df.groupby(
        [END_USE_CAT_NAMES, carrier_names_df['hh_carrier_name'].dropna().to_dict(), country_codes],
        level=['cat_code', 'carrier_code', 'country_code']
    ).sum()


    hh_end_use_df.index = hh_end_use_df.index.rename(
        ['end_use', 'carrier_name'], level=['cat_code', 'carrier_code']
    )

    # Add Swiss data
    ch_hh_end_use_df = ch_hh_consumption(ch_end_use)
    hh_end_use_df = hh_end_use_df.append(ch_hh_end_use_df, sort=True)

    # Clean up data
    hh_end_use_df = (
        hh_end_use_df
        .sort_index()
        .where(hh_end_use_df > 0)
        .dropna(how='all')
        .assign(cat_name='household')
        .set_index('cat_name', append=True)
    )

    return hh_end_use_df


def update_renewable_energy_consumption(df):
    """
    Some household consumption data has a higher overall renewables energy consumption (RA000)
    than the sum of renewable energy carriers. Here we scale all renewable energy carriers
    evenly, to match the total of RA000.
    """
    def _get_rows_to_update(df):
        renewables = (
            df
            .stack()
            .unstack("carrier_code")
            .filter(regex="^R")
            .where(lambda x: x > 0)
            .dropna(how='all')
        )
        renewables_carriers = renewables.drop('RA000', axis=1).sum(axis=1, min_count=1)
        renewables_all = renewables.xs('RA000', axis=1)
        # Only update those rows where the sum of renewable energy carriers != RA000
        return renewables.loc[~np.isclose(renewables_carriers, renewables_all)], renewables_carriers

    to_update, renewables_carriers = _get_rows_to_update(df)
    # Some rows have no data other than RA000, so we need to assign that data to one of the
    # renewable energy carriers. We choose biofuels here (R5110-5150_W6000RI)
    completely_missing = renewables_carriers[renewables_carriers.isna()].index.intersection(to_update.index)
    to_update.loc[completely_missing, 'R5110-5150_W6000RI'] = (
        to_update.loc[completely_missing, 'R5110-5150_W6000RI']
        .fillna(to_update.loc[completely_missing, 'RA000'])
    )

    # Now we scale all renewable energy carriers to match RA000
    mismatch = to_update.xs('RA000', axis=1).div(to_update.drop('RA000', axis=1).sum(axis=1))
    updated = to_update.drop('RA000', axis=1).mul(mismatch, axis=0)
    assert np.allclose(updated.sum(axis=1), to_update.xs('RA000', axis=1))

    updated_reordered = updated.stack().unstack('year').reorder_levels(df.index.names)
    # Add new rows
    df = df.append(updated_reordered.loc[updated_reordered.index.difference(df.index)])
    # Update existing rows
    df.update(updated_reordered)
    # Ensure everything has been updated as expected
    assert _get_rows_to_update(df)[0].empty


def ch_hh_consumption(ch_end_use):
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """
    ch_hh_end_use_df_sh = get_ch_sheet(
        ch_end_use, 'Tabelle 18', skipfooter=8, translation=CH_ENERGY_CARRIER_TRANSLATION
    )
    ch_hh_end_use_df_hw = get_ch_sheet(
        ch_end_use, 'Tabelle 20', skipfooter=5, translation=CH_ENERGY_CARRIER_TRANSLATION
    )
    # Quirk of the excel is that there is no space in this sheet name
    ch_hh_end_use_df_c = get_ch_sheet(
        ch_end_use, 'Tabelle21', skipfooter=4, translation=CH_ENERGY_CARRIER_TRANSLATION
    )

    ch_hh_end_use_df = (
        pd.concat(
            [ch_hh_end_use_df_sh, ch_hh_end_use_df_hw, ch_hh_end_use_df_c],
            keys=('space_heat', 'water_heat', 'cooking'),
            names=['cat_name', 'carrier_name']
        )
        .assign(country_code='CHE')
        .set_index('country_code', append=True)
    )

    # Columns are years
    ch_hh_end_use_df.columns = ch_hh_end_use_df.columns.astype(int).rename('year')

    return ch_hh_end_use_df.apply(utils.pj_to_twh)  # PJ -> TWh


def map_jrc_to_eurostat(energy_balance, jrc_end_use_df):
    jrc_end_use_df = (
        jrc_end_use_df
        .xs('ktoe', level='unit')
        .rename(utils.get_alpha3, level='country_code')
        .apply(utils.ktoe_to_twh)  # kTOE -> TWh
    )
    jrc_end_use_percent = (
        jrc_end_use_df
        .div(jrc_end_use_df.unstack('end_use').sum(axis=1))
        .unstack('year')
        .dropna(how='all')
    )
    # Fill newer, missing years
    missing_years = (
        energy_balance.columns
        .difference(jrc_end_use_df.index.get_level_values('year').unique())
    )
    for yr in missing_years:
        jrc_end_use_percent[int(yr)] = jrc_end_use_percent.mean(axis=1)

    # Fill missing countries with neighbouring countries
    jrc_end_use_percent_by_country = jrc_end_use_percent.stack('year').unstack('country_code')
    balkan_countries = jrc_end_use_percent_by_country[['BGR', 'HRV', 'HUN', 'ROU', 'GRC']].mean(axis=1)
    nordic_countries = jrc_end_use_percent_by_country[['SWE', 'FIN', 'DNK']].mean(axis=1)
    jrc_end_use_percent_by_country = jrc_end_use_percent_by_country.assign(
        ALB=balkan_countries,
        BIH=balkan_countries,
        MNE=balkan_countries,
        MKD=balkan_countries,
        SRB=balkan_countries,
        NOR=nordic_countries,
        ISL=nordic_countries
    )

    mapped_end_uses = (
        jrc_end_use_percent_by_country.stack()
        .align(energy_balance.stack())[1]
        .mul(jrc_end_use_percent_by_country.stack())
        .dropna()
    )

    return mapped_end_uses


def get_commercial_energy_consumption(
    energy_balance, ch_end_use, jrc_end_use,
    country_codes, annual_consumption
):
    """
    Use JRC IDEES service sector consumption to estimate consumption for
    space heating and hot water in the commercial sector across all countries.
    Add Swiss data from Swiss govt. stats.
    """

    # 'fuel' is actually just generic non-electric energy, which distribute based on household data
    ch_con_fuel = ch_non_hh_consumption(ch_end_use, 'Tabelle 25', annual_consumption)
    ch_con_elec = ch_non_hh_consumption(ch_end_use, 'Tabelle26', 'electricity')

    jrc_end_use_df = utils.read_tdf(jrc_end_use).xs('consumption', level='energy')

    # Map JRC end uses to annual commercial demand
    mapped_end_uses = map_jrc_to_eurostat(energy_balance, jrc_end_use_df)

    # Add Swiss data and ambient heat from heat pumps
    mapped_end_uses = (
        mapped_end_uses
        .append(ch_con_fuel.rename({'process_heat': 'cooking'}).reorder_levels(mapped_end_uses.index.names))
        .append(ch_con_elec.rename({'process_heat': 'cooking'}).reorder_levels(mapped_end_uses.index.names))
        .append(
            energy_balance  # JRC data only refers to heat pumps for heating in space heating
            .loc[['ambient_heat']]
            .assign(end_use='space_heat')
            .set_index('end_use', append=True)
            .stack()
            .reorder_levels(mapped_end_uses.index.names)
        )
    )
    mapped_end_uses.index = mapped_end_uses.index.remove_unused_levels()

    annual_consumption = annual_consumption.append(
        mapped_end_uses
        .where(mapped_end_uses > 0)
        .dropna()
        .unstack('year')
        [annual_consumption.columns]
        .assign(cat_name='commercial')
        .set_index('cat_name', append=True)
        .reorder_levels(annual_consumption.index.names)
    ).sort_index()

    return annual_consumption


def ch_non_hh_consumption(ch_end_use, sheet_name, carrier):
    """
    Switzerland data isn't in Eurostat, so we get it from their govt. stats directly
    """

    ch_con = get_ch_sheet(
        ch_end_use, sheet_name, skipfooter=4, translation=CH_HH_END_USE_TRANSLATION
    )
    # this is actually just generic non-electric energy,
    # which we assign to fuels using household ratios
    if isinstance(carrier, str):
        ch_con = ch_con.assign(carrier_name=carrier)
    else:
        hh_con = carrier.xs(('CHE', 'household'), level=('country_code', 'cat_name'))
        hh_ratios = hh_con.div(hh_con.drop('electricity', level='carrier_name').sum(level=['end_use']), axis=0)
        ch_con_disaggregated = hh_ratios.mul(ch_con, level='end_use', axis=0).dropna(how='all')

        assert np.allclose(ch_con_disaggregated.sum(level='end_use'), ch_con)
        ch_con = ch_con_disaggregated.reset_index('carrier_name')
    ch_con = (
        ch_con
        .assign(country_code='CHE')
        .set_index(['country_code', 'carrier_name'], append=True)
        .stack()
        .rename_axis(index=['end_use', 'country_code', 'carrier_name', 'year'])
        .apply(utils.pj_to_twh)  # PJ -> TWh
    )

    return ch_con


def hardcoded_country_cleanup(annual_consumption, energy_balance_df):
    # Hardcoding to ignore Iceland (known >95% geothermal), which isn't in the
    # database of household energy consumption.
    # See https://nea.is/the-national-energy-authority/energy-data/data-repository/ OS-2019-T007-01
    annual_consumption.loc[idx[:, :, 'ISL', 'household'], :] = 0

    # Montenegro data isn't in the database, so we combine estimates from the World bank
    # with average end use contributions of neighbouring countries
    # biomass fuel = 69%, electricity = 28%, oil/solid fossil fuel = 1-2%
    # See The World Bank Montenegro Second Energy Efficiency Project (P165509)
    MNE_energy_balance = energy_balance_df.loc[idx[:, 'MNE'], :]
    MNE_heat_electricity_consumption = MNE_energy_balance.loc['biofuel'] * 0.28 / 0.69
    MNE_energy_balance.loc['electricity'].update(MNE_heat_electricity_consumption)
    neighbours = ['SRB', 'HRV', 'ALB', 'BIH']
    neighbour_consumption = annual_consumption.loc[idx[END_USE_CAT_NAMES.values(), :, neighbours, ['household']], :]
    end_use_contributions = (
        neighbour_consumption
        .sum(level=['end_use', 'country_code', 'cat_name'])
        .div(neighbour_consumption.sum(level='country_code'))
        .mean(level=['end_use', 'cat_name'])
    )
    MNE_end_use = (
        end_use_contributions
        .stack().unstack(['end_use', 'cat_name'])
        .mul(MNE_energy_balance.stack(), axis=0)
        .stack([0, 1]).unstack('year')
        .reorder_levels(annual_consumption.index.names)
    )
    MNE_end_use = MNE_end_use.where(MNE_end_use > 0).dropna(how='all')
    annual_consumption = (annual_consumption.append(MNE_end_use, sort=True).sort_index())

    return annual_consumption


def fill_data_gaps(end_use_df, energy_balance_dfs, fill):

    end_use_df = end_use_df.where(end_use_df > 0)

    # Fill the household sector's end-use data based on total sectoral consumption
    hh_country_energy_balance = (
        energy_balance_dfs['hh']
        .sum(level='country_code')
        .stack()
        .where(lambda x: x > 0)
    )

    end_use_df.loc[:, idx[:, 'household']] = (
        end_use_df.loc[:, idx[:, 'household']].fillna(
            end_use_df.loc[:, idx[:, 'household']]
            .div(hh_country_energy_balance, axis=0)
            .mean(level='country_code')
            .mul(hh_country_energy_balance, level='country_code', axis=0)
        )
    )

    end_use_df = end_use_df.where(end_use_df > 0)

    # For all remaining gaps, fill with mean/first/last/max/min/whatever data for the
    # country, based on the string given by 'fill'
    end_use_df = end_use_df.fillna(
        end_use_df.groupby('country_code').agg(fill)
    )

    return end_use_df


def get_annual_electricity_consumption(annual_consumption, energy_balance_dfs):
    """
    Get annual energy consumption coming from electricity, to remove later from the
    electricity demand profile.
    Gaps in electricity consumption are filled before saving, based on annual energy
    consumption
    """
    electricity_consumption = (
        annual_consumption
        .drop('end_use_electricity')
        .loc[idx[:, ['electricity', 'direct_electric', 'heat_pump'], :, :], :]
        .sum(level=['end_use', 'country_code', 'cat_name'])
        .stack().unstack(['end_use', 'cat_name'])
    )

    electricity_consumption = fill_data_gaps(
        electricity_consumption, energy_balance_dfs, fill='first'
    )
    electricity_consumption = (
        electricity_consumption
        .rename(lambda x: f'{x}_bau_electricity', level='end_use', axis=1)
    )
    return electricity_consumption


def get_national_heat_demand(annual_consumption, energy_balance_dfs, heat_tech_params):
    """
    Get consumption of heat techs and convert to demand for heating
    """
    def efficiencies(end_use):
        return pd.Series({
            'biogas': heat_tech_params[f'{end_use}_techs'].get('gas_eff', 1),
            'biofuel': heat_tech_params[f'{end_use}_techs'].get('biofuel_eff', 1),
            'solid_fossil': heat_tech_params[f'{end_use}_techs'].get('solid_fossil_eff', 1),
            'natural_gas': heat_tech_params[f'{end_use}_techs'].get('gas_eff', 1),
            'manufactured_gas': heat_tech_params[f'{end_use}_techs'].get('gas_eff', 1),
            'gas': heat_tech_params[f'{end_use}_techs'].get('gas_eff', 1),
            'oil': heat_tech_params[f'{end_use}_techs'].get('oil_eff', 1),
            'solar_thermal': heat_tech_params[f'{end_use}_techs'].get('solar_thermal_eff', 1),
            'renewable_heat': heat_tech_params[f'{end_use}_techs'].get('solar_thermal_eff', 1),
            'electricity': heat_tech_params[f'{end_use}_techs'].get('electricity_eff', 1),
            'direct_electric': 1,  # don't need to deal with heat pump COP if direct electric is 100% efficient
            'heat': 1,
            # heat demand met by heat pumps = heat pump electricity + ambient heat
            'heat_pump': 1,
            'ambient_heat': 1
        })

    # Convert end use consumption data to demand
    demands = []
    for end_use in ['space_heat', 'water_heat', 'cooking']:
        _demand = (
            annual_consumption.loc[[end_use]]
            .mul(efficiencies(end_use), level='carrier_name', axis=0)
            .sum(level=['end_use', 'country_code', 'cat_name'])
        )
        # sense check: demand is at least the minimum efficiency * consumption
        assert (_demand >= (
            annual_consumption.loc[[end_use]]
            .mul(efficiencies(end_use).min())
            .sum(level=['end_use', 'country_code', 'cat_name'])
        )).all().all()
        demands.append(_demand)

    demand = pd.concat(demands).stack().unstack(['end_use', 'cat_name'])

    # Fill gaps in demand using energy balance data
    # This is done by finding the average contribution of cooking, space heating
    # and water heating to HH demand (from energy balances),
    # then applying to years in which HH sub class data is not available.
    demand = fill_data_gaps(
        demand,
        energy_balance_dfs,
        fill='mean'
    )

    return demand


def get_ch_sheet(path_to_excel, sheet, skipfooter, translation=None):
    _df = (
        pd.read_excel(
            path_to_excel, sheet_name=sheet,
            skiprows=9, skipfooter=skipfooter, index_col=1
        )
        .drop(['Unnamed: 0', 'Δ ’00 – ’18'], axis=1)
    )
    _df.index = _df.index.str.strip()
    _df.columns = _df.columns.astype(int)

    if translation is not None:
        return _df.groupby(translation).sum()
    else:
        return _df


if __name__ == "__main__":
    get_heat_demand(
        path_to_hh_end_use=snakemake.input.hh_end_use,
        path_to_ch_end_use=snakemake.input.ch_end_use,
        path_to_energy_balance=snakemake.input.energy_balance,
        path_to_commercial_demand=snakemake.input.commercial_demand,
        path_to_carrier_names=snakemake.input.carrier_names,
        countries=snakemake.params.countries,
        heat_tech_params=snakemake.params.heat_tech_params,
        path_to_electricity_consumption=snakemake.output.electricity,
        path_to_output=snakemake.output.demand,
    )