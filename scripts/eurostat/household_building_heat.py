import pandas as pd
import numpy as np

from eurocalliopelib import utils

YEARS = slice(2000, 2018)


def household_end_use_energy(
    path_to_household_end_use_energy_balance: str, path_to_output: str
):
    """
    Extract household end use energy consumption from Eurostat
    """

    # Index name in TSV file is 'nrg_bal,siec,unit,geo\time'
    hh_end_use_df = utils.read_eurostat_tsv(path_to_household_end_use_energy_balance)

    index_names = {
        "nrg_bal": "cat_code",
        "siec": "carrier_code",
        "geo": "country_code",
        "time": "year"
    }
    # Just keep one unit's worth of data (we choose TJ)
    hh_end_use_dataarray = (
        hh_end_use_df
        .xs('TJ', level='unit')
        .apply(utils.tj_to_twh)  # TJ -> TWh
        .astype(float)
        .dropna(how='all')
        .stack()
        .to_xarray()
        .rename(index_names)
        .sortby("year")
    )

    # clean up renewables info
    hh_end_use_dataarray = update_renewable_energy_consumption(hh_end_use_dataarray)

    # Add missing renewables data to

    country_code_mapping = utils.convert_valid_countries(hh_end_use_dataarray.country_code.values)
    hh_end_use_dataarray = utils.rename_and_groupby(hh_end_use_dataarray, country_code_mapping, dim_name="country_code")

    # Clean up data
    hh_end_use_dataarray = hh_end_use_dataarray.sel(year=YEARS) # limit year of data to 2018 as the most recent

    hh_end_use_dataarray.assign_attrs(unit="twh").to_netcdf(path_to_output)


def update_renewable_energy_consumption(energy_consumption_da):
    """
    Some household consumption data has a higher overall renewables energy consumption (RA000)
    than the sum of renewable energy carriers. Here we scale all renewable energy carriers
    evenly, to match the total of RA000.
    """
    energy_consumption_da_to_check = energy_consumption_da.copy()

    renewables_sub_carriers, renewables_sub_carriers_sum, renewables_all = get_renewables_carrier_data(energy_consumption_da)

    # Some data points have no renewables data other than the sum (RA000),
    # so we need to assign that data to one of the renewable energy sub-carriers.
    # ASSUME: if no data in any renewables sub-carriers, energy is from biofuels (R5110-5150_W6000RI)
    energy_consumption_da.loc[{"carrier_code": 'R5110-5150_W6000RI'}] = (
        energy_consumption_da
        .where(energy_consumption_da > 0)
        .sel(carrier_code='R5110-5150_W6000RI')
        .fillna(renewables_all.where(renewables_sub_carriers_sum.isnull() & renewables_all.notnull()))
    )

    coords_to_update = renewables_sub_carriers.to_series().dropna().to_xarray().coords
    required_renewables_scaling = renewables_all / renewables_sub_carriers_sum
    energy_consumption_da.loc[coords_to_update] = (
        energy_consumption_da.sel(coords_to_update) * required_renewables_scaling.fillna(1)
    )

    _, new_renewables_sub_carriers_sum, new_renewables_all = get_renewables_carrier_data(energy_consumption_da)
    _, _, orig_renewables_all = get_renewables_carrier_data(energy_consumption_da_to_check)

    # Test that the sum of renewable energy carriers equals the total defined by RA000
    assert np.allclose(new_renewables_sub_carriers_sum, new_renewables_all, equal_nan=True)
    assert np.allclose(orig_renewables_all, new_renewables_all, equal_nan=True)
    assert np.allclose(
        energy_consumption_da_to_check.where(~energy_consumption_da.carrier_code.str.startswith("R")),
        energy_consumption_da.where(~energy_consumption_da.carrier_code.str.startswith("R")),
        equal_nan=True
    )

    return energy_consumption_da


def get_renewables_carrier_data(energy_consumption_da):
    renewables_carriers = (
        energy_consumption_da
        .where(energy_consumption_da.carrier_code.str.startswith("R") & (energy_consumption_da > 0))
    )
    renewables_sub_carriers = renewables_carriers.drop_sel(carrier_code="RA000")
    renewables_sub_carriers_sum = renewables_sub_carriers.sum("carrier_code", min_count=1)
    renewables_all = renewables_carriers.sel(carrier_code="RA000")

    return renewables_sub_carriers, renewables_sub_carriers_sum, renewables_all


if __name__ == "__main__":
    household_end_use_energy(
        path_to_household_end_use_energy_balance=snakemake.input.household_end_use_energy_balance,
        path_to_output=snakemake.output[0]
    )
