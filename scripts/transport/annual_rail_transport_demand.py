import pandas as pd

from eurocalliopelib import utils

from annual_road_transport_demand import get_all_distance_efficiency

idx = pd.IndexSlice


def get_rail_transport_demand(
    energy_balances_path, jrc_rail_energy_path, jrc_rail_distance_path,
    rail_energy_out_path, rail_bau_electricity_out_path
):
    energy_balances = utils.read_tdf(energy_balances_path)
    rail_energy_df = utils.read_tdf(jrc_rail_energy_path)
    rail_distance_df = utils.read_tdf(jrc_rail_distance_path)

    _total_rail_distance, rail_efficiency, rail_bau_consumption = get_all_distance_efficiency(
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
    )
    rail_energy_twh = utils.add_idx_level(rail_energy, unit="twh")

    rail_electricity_bau = (
        rail_bau_consumption
        .sum(level=['carrier', 'section', 'country_code', 'year', 'unit'])
        .xs('electricity')
    )
    # High speed and metro are all electrified, so bau_electricity == eurocalliope_electricity
    assert abs(
        (rail_energy_twh - rail_bau_consumption.xs('electricity', level='carrier'))
        [['High speed passenger trains', 'Metro and tram, urban light rail']]
    ).max() < 1e-10

    rail_energy_twh.to_csv(rail_energy_out_path)
    rail_electricity_bau.to_csv(rail_bau_electricity_out_path)


if __name__ == "__main__":
    get_rail_transport_demand(
        energy_balances_path=snakemake.input.energy_balances,
        jrc_rail_energy_path=snakemake.input.jrc_rail_energy,
        jrc_rail_distance_path=snakemake.input.jrc_rail_distance,
        rail_energy_out_path=snakemake.output.rail_energy,
        rail_bau_electricity_out_path=snakemake.output.rail_bau_electricity,
    )
