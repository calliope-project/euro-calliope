import pandas as pd

from eurocalliopelib import utils

idx = pd.IndexSlice


def get_air_transport_demand(
    energy_balances_path, air_energy_out_path,
):
    energy_balances = utils.read_tdf(energy_balances_path)

    # Add air transport energy demand from 'not elsewhere specified' (military)
    # ASSUME: 'not elsewhere specified' oil use goes predominantly to 'road' transport,
    # except kerosene/jet-fuel, which goes to aviation
    other_transport_aviation = (  # all kerosene from the military destined for aviation
        energy_balances
        .loc[idx['FC_OTH_NSP_E', ['O4651', 'O4653', 'O4661XR5230B', 'O4669'], :, :, :]]
        .sum(level=['country', 'year'])
    )

    air_energy = (
        energy_balances
        .loc[idx[['FC_TRA_DAVI_E', 'INTAVI'], 'O4000XBIO', :, :, :]]
        .sum(level=['country', 'year'])
        .add(other_transport_aviation, fill_value=0)
    )
    air_energy_twh = utils.add_idx_level(air_energy.apply(utils.tj_to_twh), unit="twh")

    air_energy_twh.to_csv(air_energy_out_path)


if __name__ == "__main__":
    get_air_transport_demand(
        energy_balances_path=snakemake.input.energy_balances,
        air_energy_out_path=snakemake.output.air_energy,
    )
