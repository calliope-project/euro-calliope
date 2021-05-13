import pandas as pd

from eurocalliopelib import utils

CARRIERS = {
    'O4652XR5210B': 'petrol',
    'R5210B': 'biofuels',
    'O4671XR5220B': 'diesel',
    'R5220B': 'biofuels',
    'O4630': 'lpg',
    'G3000': 'natural_gas',
    'E7000': 'electricity'
}

idx = pd.IndexSlice


def get_rail_transport_demand(
    energy_balances_path, jrc_rail_energy_path, jrc_rail_distance_path, 
    rail_energy_out_path, rail_bau_electricity_out_path
):
    energy_balances = utils.read_tdf(energy_balances_path)
    rail_energy_df = utils.read_tdf(jrc_rail_energy_path)
    rail_distance_df = utils.read_tdf(jrc_rail_distance_path)

    total_rail_distance, rail_efficiency, rail_bau_consumption = get_all_distance_efficiency(
        energy_balances, 'FC_TRA_RAIL_E', rail_energy_df, rail_distance_df, 'carrier'
    )
    # Some cleanup that's specific to rail data
    rail_energy = (
        # historical energy consumption -> historical distance
        rail_bau_consumption.droplevel('unit')
        .mul(rail_efficiency.droplevel('unit'))    # efficiency = distance/energy
        .sum(level=['vehicle_type', 'section', 'country_code', 'year'])
        # historical distance -> total electricity demand
        .div(rail_efficiency.xs('electricity', level='carrier').droplevel('unit'))
        .sum(level=['vehicle_type', 'country_code', 'year'])
        .to_frame('twh')
        .rename_axis(columns='unit')
        .stack()
    )
    rail_electricity_bau = (
        rail_bau_consumption
        .sum(level=['carrier', 'section', 'country_code', 'year', 'unit'])
        .xs('electricity')
    )
    # High speed and metro are all electrified, so bau_electricity == eurocalliope_electricity
    assert abs(
        (rail_energy - rail_bau_consumption.xs('electricity', level='carrier'))
        [['High speed passenger trains', 'Metro and tram, urban light rail']]
    ).max() < 1e-10

    rail_energy.to_csv(rail_energy_out_path)
    rail_electricity_bau.to_csv(rail_bau_electricity_out_path)


def get_all_distance_efficiency(
    energy_balances, cat_name, energy_df, distance_df, unique_dim, other_transport_road=0
):

    transport_energy_balance = (
        energy_balances
        .xs(cat_name)
        .unstack('carrier_code')
        .groupby(CARRIERS, axis=1).sum(min_count=1)
        .rename_axis(columns='carrier')
        .droplevel('unit')
        .stack()
        .add(other_transport_road, fill_value=0)
        .apply(utils.tj_to_twh)
        .rename_axis(index=['country_code', 'year', 'carrier'])
    )

    # contribution of each transport mode to carrier consumption from JRC_IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    carrier_contribution = fill_missing_countries_and_years(
        energy_df.div(energy_df.sum(level=['carrier', 'country_code', 'year']))
    )
    # Energy consumption per transport mode by mapping transport mode
    # carrier contributions to total carrier consumption
    transport_energy_per_mode = pd.merge(
        transport_energy_balance.to_frame('eurostat'),
        carrier_contribution.to_frame('jrc'),
        left_index=True,
        right_on=transport_energy_balance.index.names
    )
    transport_energy_per_mode = (
        transport_energy_per_mode['eurostat'].mul(transport_energy_per_mode['jrc'])
    )
    # Distance per unit energy consumed per transport mode according to JRC IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    transport_efficiency = fill_missing_countries_and_years(
        distance_df
        .where(distance_df > 0)
        .div(
            energy_df
            .where(energy_df > 0)
            .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        )
    )
    # Distance travelled per transport mode, including years 2016-2018,
    # based on JRC IDEES transport efficiency (2015 data for 2016-2018)
    transport_distance_all_years = (
        transport_energy_per_mode
        .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        .mul(transport_efficiency)
    )
    # Use the distance traveled based on Eurostat to fill in blanks in JRC data
    aligned_dfs = (
        transport_distance_all_years
        .reorder_levels(distance_df.index.names)
        .align(distance_df)
    )
    total_transport_distance = (
        aligned_dfs[1]
        .fillna(aligned_dfs[0])
        .sum(level=['vehicle_type', unique_dim, 'country_code', 'unit', 'year'])
    )

    return total_transport_distance, transport_efficiency, transport_energy_per_mode

def fill_missing_countries_and_years(jrc_data):
    jrc_data = jrc_data.unstack('country_code')
    balkan_countries = jrc_data[['BG', 'HR', 'HU', 'RO', 'EL']].mean(axis=1)
    nordic_countries = jrc_data[['SE', 'DK']].mean(axis=1)
    ch_neighbours = jrc_data[['DE', 'AT', 'FR', 'IT']].mean(axis=1)
    jrc_data = jrc_data.assign(
        AL=balkan_countries,
        BA=balkan_countries,
        ME=balkan_countries,
        MK=balkan_countries,
        RS=balkan_countries,
        NO=nordic_countries,
        IS=nordic_countries,
        CH=ch_neighbours,
    ).stack().unstack('year')
    jrc_data = jrc_data.assign(
        **{str(i): jrc_data[2015] for i in range(2016, 2019)}
    )
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()
    

if __name__ == "__main__":
    get_rail_transport_demand(
        energy_balances_path=snakemake.input.energy_balances,
        jrc_rail_energy_path=snakemake.input.jrc_rail_energy,
        jrc_rail_distance_path=snakemake.input.jrc_rail_distance,
        rail_energy_out_path=snakemake.output.rail_energy,
        rail_bau_electricity_out_path=snakemake.output.rail_bau_electricity,
    )