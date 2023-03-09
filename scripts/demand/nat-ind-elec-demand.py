import pandas as pd
import pycountry

from eurocalliopelib import utils

ENGINE = "pyxlsb"


def extract_ind_elec_demand(path_energy_balances_foldername: str,
                            path_energy_balances_filename: str,
                            countries: list,
                            year: str,
                            output_path: str) -> None:
    """
    Extracts electricity demand of industrial sectors for each country, merges to dataframe, saves as csv.

    Parameters:
        paths_energy_balances (list):
            .xlsb files with EUROSTAT energy balances for each country
        countries_alpha_3 (list):
            List of country codes in local scope in alpha_3 format
        year (str):
            Determines year of industrial demand used, i.e., datasheet that is exported.
        output_path (str):
            Location for output csv file.
    """

    df_ind_elec_demand = pd.DataFrame()
    for country in countries:
        country_code_alpha_3 = pycountry.countries.lookup(country).alpha_3
        energy_balance_path = path_energy_balances_foldername + utils.iso3_to_eu_country_code(
            country_code_alpha_3) + path_energy_balances_filename
        df_energy_balance = pd.read_excel(energy_balance_path, engine=ENGINE, sheet_name=year,
                                          header=4, usecols="H:CB", index_col=0, skipfooter=13)
        df_ind_elec_demand[country_code_alpha_3] = df_energy_balance.loc["FC_IND_IS_E":"FC_IND_NSP_E", "Electricity"]

    df_ind_elec_demand.replace("Z", 0, inplace=True) # dataset uses Z for missing data, replace by 0
    df_ind_elec_demand.to_csv(output_path)


if __name__ == "__main__":
    extract_ind_elec_demand(
        path_energy_balances_foldername=snakemake.params.path_energy_balances_foldername,
        path_energy_balances_filename=snakemake.params.path_energy_balances_filename,
        countries=snakemake.params.countries,
        year=snakemake.params.year,
        output_path=snakemake.output[0]
    )
