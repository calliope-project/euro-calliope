import pandas as pd

from eurocalliopelib import utils

idx = pd.IndexSlice


def get_marine_transport_demand(
    energy_balances_path, marine_energy_out_path
):
    energy_balances = utils.read_tdf(energy_balances_path)
    
    # Add marine transport energy demand from agriculture
    # assumption: fisheries' oil use goes to marine fuel demand;
    
    other_transport_marine = (  # i.e. all fisheries' oil use
        energy_balances
        .loc[idx[['FC_OTH_FISH_E'], 'O4000XBIO', :, :, :]]
        .sum(level=['country', 'year'])
    )

    marine_energy = (
        energy_balances
        .loc[['INTMARB', 'FC_TRA_DNAVI_E']]
        .unstack('carrier_code')
        ['O4000XBIO']
        .sum(level=['country', 'year'])
        .add(other_transport_marine, fill_value=0)
        .apply(utils.tj_to_twh)
        .to_frame('twh')
        .rename_axis(columns='unit')
        .stack()
    )

    marine_energy.to_csv(marine_energy_out_path)


if __name__ == "__main__":
    get_marine_transport_demand(
        energy_balances_path=snakemake.input.energy_balances,
        marine_energy_out_path=snakemake.output.marine_energy,
    )