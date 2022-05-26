import pandas as pd
import xarray as xr

from eurocalliopelib import utils

END_USE_CAT_NAMES = {
    'FC_OTH_HH_E_CK': 'cooking',
    'FC_OTH_HH_E_SH': 'space_heat',
    'FC_OTH_HH_E_WH': 'water_heat'
}
MNE_NEIGHBOURS = ['SRB', 'HRV', 'ALB', 'BIH']


def household_heat_energy_supply(
    path_to_ch_building_heat_energy_supply: str,
    path_to_annual_energy_balances: str,
    path_to_eurostat_household_building_heat: str,
    path_to_carrier_names: str,
    countries: list,
    path_to_output: str
):

    carrier_name_mapping = pd.read_csv(path_to_carrier_names, index_col=0, header=0)

    household_energy_balance = xr.open_dataarray(path_to_annual_energy_balances)

    # Get annual energy balance data for households
    household_end_use_energy_balance = xr.open_dataarray(
        path_to_eurostat_household_building_heat
    )

    # rename index labels to be more readable
    household_end_use_energy_balance = utils.rename_and_groupby(
        household_end_use_energy_balance, END_USE_CAT_NAMES,
        dim_name="cat_code", new_dim_name="end_use"
    )
    household_end_use_energy_balance = utils.rename_and_groupby(
        household_end_use_energy_balance,
        carrier_name_mapping['household_carrier_name'].dropna().to_dict(),
        dim_name="carrier_code", new_dim_name="carrier_name"
    )

    ch_household_end_use_energy_balance = xr.open_dataarray(path_to_ch_building_heat_energy_supply)

    household_end_use_energy_balance_incl_ch = utils.merge_da([
        household_end_use_energy_balance.expand_dims(cat_name=["household"]),
        ch_household_end_use_energy_balance.sel(cat_name=["household"])
    ])

    # fill missing data for countries with no household end-use energy supply data
    household_end_use_energy_balance_incl_ch = fill_missing_country_data_gaps(
        household_end_use_energy_balance_incl_ch, household_energy_balance
    )

    country_codes = [utils.convert_country_code(country) for country in countries]
    household_end_use_energy_balance_incl_ch = (
        household_end_use_energy_balance_incl_ch.sel(country_code=country_codes)
    )

    household_end_use_energy_balance_incl_ch.to_netcdf(path_to_output)


def fill_missing_country_data_gaps(household_end_use_energy_balance, household_annual_energy_balance):
    # ASSUME: Ignoring building heat demand in Iceland (known >95% geothermal),
    # which isn't in the database of household end-use energy balances.
    # See https://nea.is/the-national-energy-authority/energy-data/data-repository/ OS-2019-T007-01

    household_end_use_energy_balance = utils.merge_da([
        household_end_use_energy_balance,
        household_end_use_energy_balance.sum("country_code").clip(max=0).expand_dims(country_code=["ISL"])
    ])

    # ASSUME: MNE fuel consumption for building heat are distributed as:
    # biomass fuel = 69%, electricity = 28%, oil/solid fossil fuel = 1-2% (The World Bank Montenegro Second Energy Efficiency Project (P165509).

    # ASSUME: MNE distribution of building heat energy supply to end uses has the same
    # relative distribution as in neighbouring countries.

    # ASSUME: all biofuel consumption in MNE is for building heat
    MNE_energy_balance = household_annual_energy_balance.sel(country_code="MNE")

    MNE_heat_electricity_consumption = MNE_energy_balance.sel(carrier_name="biofuel") * 0.28 / 0.69
    MNE_energy_balance.loc[{"carrier_name": "electricity"}] = MNE_heat_electricity_consumption

    neighbour_consumption = household_end_use_energy_balance.sel(country_code=MNE_NEIGHBOURS)
    neighbour_end_use_contribution_fractions = (
        (neighbour_consumption.sum("carrier_name") / neighbour_consumption.sum("country_code"))
        .mean("country_code")
    )
    MNE_end_use_energy_balance = MNE_energy_balance * neighbour_end_use_contribution_fractions

    return utils.merge_da([household_end_use_energy_balance, MNE_end_use_energy_balance.expand_dims(country_code=["MNE"])])


if __name__ == "__main__":
    household_heat_energy_supply(
        path_to_ch_building_heat_energy_supply=snakemake.input.ch_building_heat_energy_supply,
        path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
        path_to_eurostat_household_building_heat=snakemake.input.eurostat_household_building_heat,
        path_to_carrier_names=snakemake.params.carrier_names,
        countries=snakemake.params.countries,
        path_to_output=snakemake.output[0],
    )
