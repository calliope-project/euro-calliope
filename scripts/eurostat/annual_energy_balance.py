import pandas as pd

from eurocalliopelib import utils

YEARS = range(2000, 2018)


def generate_annual_energy_balance_nc(
    path_to_input,
    path_to_cat_names,
    path_to_carrier_names,
    path_to_result,
):
    """
    Open a TSV file and reprocess it into a xarray dataset, including long names for
    Eurostat codes.
    Switzerland is not included in Eurostat, so we splice in data from their govt.
    statistics.
    """
    # Names for each consumption category/sub-category and carriers have been prepared manually
    cat_names = pd.read_csv(path_to_cat_names, header=0, index_col=0)
    carrier_names = pd.read_csv(path_to_carrier_names, header=0, index_col=0)

    index_names = {
        "nrg_bal": "cat_code",
        "siec": "carrier_code",
        "geo": "country_code",
        "time": "year"
    }
    da = utils.read_eurostat_tsv(path_to_input).stack().to_xarray().rename(index_names)

    da = da.sel(
        cat_code=cat_names.index.intersection(da.cat_code),
        carrier_code=carrier_names.index.intersection(da.carrier_code),
        unit="TJ"
    )

    country_code_mapping = utils.convert_valid_countries(da.country_code.values)
    da = utils.rename_and_groupby(da, country_code_mapping, dim_name="country_code")

    da = utils.tj_to_twh(da).drop_vars("unit").assign_attrs({"unit": "twh"})

    da = da.sel(year=YEARS).rename("annual_energy_balances")

    da.to_netcdf(path_to_result)


if __name__ == "__main__":
    generate_annual_energy_balance_nc(
        path_to_input=snakemake.input.eurostat_energy_balance,
        path_to_cat_names=snakemake.params.cat_names,
        path_to_carrier_names=snakemake.params.carrier_names,
        path_to_result=snakemake.output[0],
    )
