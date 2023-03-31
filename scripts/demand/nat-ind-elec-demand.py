import pandas as pd
import pycountry

from eurocalliopelib import utils

ENGINE_xlsb = "pyxlsb"
ENGINE_xlsx = "openpyxl"


def extract_ind_elec_demand(path_energy_balances_foldername: str,
                            path_energy_balances_europe_filename: str,
                            path_energy_balances_switzerland_filename: str,
                            countries: list,
                            year: str,
                            output_path: str) -> None:
    """
    Extracts electricity demand of industrial sectors for each country, merges to dataframe, saves as csv.

    Parameters:
        paths_energy_balances_foldername (str):
            Foldername with energy balances
        path_energy_balances_europe_filename (str):
            Filenames of .xlsb files with EUROSTAT energy balances for each country but Switzerland
        path_energy_balances_europe_filename (str):
            Filenames of .xlsx files BFE energy balances for Switzerland
        countries (list):
            List of country names in scope
        year (str):
            Determines year of industrial demand used
        output_path (str):
            Location for output csv file.
    """

    df_ind_elec_demand = pd.DataFrame()
    for country in countries:
        country_code_alpha_3 = pycountry.countries.lookup(country).alpha_3

        if country != "Switzerland":
            df_ind_elec_demand = extract_from_balance_sheet_europe(country_code_alpha_3, df_ind_elec_demand, year,
                                                                   path_energy_balances_foldername,
                                                                   path_energy_balances_europe_filename)
        else:
            df_ind_elec_demand = extract_from_balance_sheet_switzerland(country_code_alpha_3, df_ind_elec_demand, year,
                                                                        path_energy_balances_foldername,
                                                                        path_energy_balances_switzerland_filename)

    df_ind_elec_demand.fillna(0, inplace=True) # NaNs occur if in a sector CH is only country with demand.
    df_ind_elec_demand.to_csv(output_path)


def extract_from_balance_sheet_europe(country_code_alpha_3: str, df_ind_elec_demand: pd.DataFrame, year: str,
                                      path_energy_balances_foldername: str,
                                      path_energy_balances_europe_filename: str) -> pd.DataFrame:
    """
    Extracts electricity demand of industrial sectors for all countries but Switzerland.

    Parameters:
        paths_energy_balances_foldername (str):
            Foldername with energy balances
        path_energy_balances_europe_filename (str):
            Filenames of .xlsb files with EUROSTAT energy balances for each country but Switzerland
        country_code_alpha_3 (str):
            Country code of country at hand, in alpha 3 format
        year (str):
            Determines year of industrial demand extracted, i.e., the sheet
        output_path (str):
            Location for output csv file.
    Returns:
        df_ind_elec_demand (pd.DataFrame):
            Returns dataframe where the industrial demand data was added for the country at hand
    """

    energy_balance_path = path_energy_balances_foldername + utils.iso3_to_eu_country_code(
        country_code_alpha_3) + path_energy_balances_europe_filename
    df_energy_balance = pd.read_excel(energy_balance_path, engine=ENGINE_xlsb, sheet_name=year,
                                      header=4, usecols="H:CB", index_col=0, skipfooter=13)
    # For now only industrial sectors, not "transport" and "others":
    df_ind_elec_demand[country_code_alpha_3] = (
        df_energy_balance.loc["FC_IND_IS_E":"FC_IND_NSP_E", "Electricity"] * 11630) # ktoe to MWh,
    df_ind_elec_demand.replace("Z", 0, inplace=True) # dataset uses Z for missing data, replace by 0

    return df_ind_elec_demand


def extract_from_balance_sheet_switzerland(country_code_alpha_3: str, df_ind_elec_demand: pd.DataFrame, year: str,
                                           path_energy_balances_foldername: str,
                                           path_energy_balances_switzerland_filename: str) -> pd.DataFrame:
    """
    Extracts electricity demand of industrial sectors for all countries but Switzerland.

    Parameters:
        paths_energy_balances_foldername (str):
            Foldername with energy balances
        path_energy_balances_switzerland_filename (str):
            Filenames of .xlsx files with BFE energy balances for Switzerland
        country_code_alpha_3 (str):
            Country code of country at hand, in alpha 3 format
        year (str):
            Determines year of industrial demand extracted, i.e., the column
        output_path (str):
            Location for output csv file.
    Returns:
        df_ind_elec_demand (pd.DataFrame):
            Returns dataframe where the industrial demand data was added for the country at hand
    """

    energy_balance_path = (path_energy_balances_foldername + utils.iso3_to_eu_country_code(country_code_alpha_3)
                           + path_energy_balances_switzerland_filename)

    df_ind_final_energy_consumption_PJ = pd.read_excel(energy_balance_path, engine=ENGINE_xlsx, sheet_name="Tabelle30",
                                                       header=4, usecols="B:X", index_col=0, skipfooter=7)
    df_ind_final_elec_consumption_PJ = pd.read_excel(energy_balance_path, engine=ENGINE_xlsx, sheet_name="Tabelle32",
                                                     header=4, usecols="B:X", index_col=0, skipfooter=7)
    df_ind_sector_shares_energy_consumption = pd.read_excel(energy_balance_path, engine=ENGINE_xlsx,
                                                            sheet_name="Tabelle33", header=4, usecols="B:H",
                                                            index_col=0, skipfooter=7)

    # Group Verwendungszwecke (Raumwärme, Prozesswärme, Mechanische Arbeit, etc.) so that same format as in
    # df_ind_sector_shares_energy_consumption_PJ
    for df in [df_ind_final_energy_consumption_PJ, df_ind_final_elec_consumption_PJ]:
        df.loc["Raumwärme & Warmwasser"] = df.loc["Raumwärme"] + df.loc["Warmwasser"]
        df.loc["Beleuchtung, HT, I&K"] = (df.loc["Beleuchtung"] + df.loc["Klima, Lüftung, HT"]
                                          + df.loc["I&K, Unterhaltung"])
        df.rename({"Antriebe, Prozesse": "Mechanische Arbeit",
                   "sonstige": "Elektrolyse, Umweltschutz und sonstige"}, inplace=True)
        df.drop(["Raumwärme", "Warmwasser", "Beleuchtung", "Klima, Lüftung, HT", "I&K, Unterhaltung"], inplace=True)
    df_ind_final_elec_consumption_PJ.drop(["Total Elektrizität"], inplace=True)
    df_ind_final_energy_consumption_PJ.drop(["Total Endenergie"], inplace=True)
    df_ind_sector_shares_energy_consumption.drop(columns=["Anteil am Energieverbrauch"], inplace=True)

    # Slicing to the relevant year
    df_ind_final_energy_consumption_PJ = df_ind_final_energy_consumption_PJ[year]
    df_ind_final_elec_consumption_PJ = df_ind_final_elec_consumption_PJ[year]

    # Computations:
    # For each Verwendungszweck (Raumwärme, Prozesswärme, Mechanische Arbeit, etc.) compute ratio of industrial swiss
    # final elec consumption / industrial swiss final energy consumption
    share_elec_of_energy = df_ind_final_elec_consumption_PJ / df_ind_final_energy_consumption_PJ
    # For each Verwendungszweck compute final energy consumption for each industrial sector using the provided shares
    df_ind_final_energy_consumption_sector_PJ = (
        df_ind_sector_shares_energy_consumption * df_ind_final_energy_consumption_PJ)
    # For each Verwendungszweck compute final electricity demand from aboves ratio and the final energy demand
    df_ind_final_elec_consumption_sector_PJ = df_ind_final_energy_consumption_sector_PJ * share_elec_of_energy
    # For each sector sum demands over all Verwendungszwecke to get swiss final ind elec demand for each sector
    df_ind_final_elec_consumption_sector_PJ[country_code_alpha_3] = df_ind_final_elec_consumption_sector_PJ.sum(axis=1)
    df_ind_final_elec_consumption_sector_PJ = df_ind_final_elec_consumption_sector_PJ[country_code_alpha_3]
    # ASSUMES that the ratio "Final electricty consumption of Switzerland's industry / Final energy consumption of
    # Switzerland's industry" euqals the ratio "Final electricty consumption of Switzerland's industrial sector j /
    # Final energy consumption of Switzerland's industrial sector j" for all industrial sectors (only needs to hold
    # within each end use category, i.e., "Verwendungszwecke" (Raumwärme, Prozesswärme, Mechanische Arbeit, etc.).

    # Group by EUROSTAT activity categories and rename
    # ASSUMES that industrial demand is negligible (=0) in sectors that are not reported by Swiss BFE (i.e.,
    # FC_IND_TE_E (transport equipment), FC_IND_WP_E (wood products), FC_IND_MQ_E (mining and quarrying))
    df_ind_final_elec_consumption_sector_PJ.loc["Glas, Keramik, Beton, Steine, Zement, Kalk, Ziegel"] = (
        df_ind_final_elec_consumption_sector_PJ.loc["Glas, Keramik, Beton, Steine"]
        + df_ind_final_elec_consumption_sector_PJ.loc["Zement, Kalk, Ziegel"])
    df_ind_final_elec_consumption_sector_PJ.loc["Maschinenbau, Fahrzeugbau, Metallerzeugnisse, Geräte"] = (
        df_ind_final_elec_consumption_sector_PJ.loc["Maschinenbau, Fahrzeugbau"]
        + df_ind_final_elec_consumption_sector_PJ.loc["Metallerzeugnisse, Geräte"])
    df_ind_final_elec_consumption_sector_PJ.loc["Übrige, Wasser, Abfall"] = (
        df_ind_final_elec_consumption_sector_PJ.loc["Übrige"]
        + df_ind_final_elec_consumption_sector_PJ.loc["Wasser, Abfall"])
    df_ind_final_elec_consumption_sector_PJ.drop(["Glas, Keramik, Beton, Steine", "Zement, Kalk, Ziegel",
                                                  "Maschinenbau, Fahrzeugbau", "Metallerzeugnisse, Geräte", "Übrige",
                                                  "Wasser, Abfall", "Total"], inplace=True)
    df_ind_final_elec_consumption_sector_PJ.rename(
        {"Bau": "FC_IND_CON_E",
         "Chemie, Pharma": "FC_IND_CPC_E",
         "Eisen, Stahl": "FC_IND_IS_E",
         "Glas, Keramik, Beton, Steine, Zement, Kalk, Ziegel": "FC_IND_NMM_E",
         "Maschinenbau, Fahrzeugbau, Metallerzeugnisse, Geräte": "FC_IND_MAC_E",
         "Nahrung, Tabak": "FC_IND_FBT_E",
         "NE-Metalle": "FC_IND_NFM_E",
         "Papier, Druck": "FC_IND_PPP_E",
         "Textilien": "FC_IND_TL_E",
         "Übrige, Wasser, Abfall": "FC_IND_NSP_E"}, inplace=True)

    df_ind_final_elec_consumption_sector_MWh = df_ind_final_elec_consumption_sector_PJ * 277777 # Convert PJ to MWh
    df_ind_elec_demand = pd.merge(df_ind_elec_demand, df_ind_final_elec_consumption_sector_MWh, left_index=True,
                                  right_index=True, how="outer")

    return df_ind_elec_demand


if __name__ == "__main__":
    extract_ind_elec_demand(
        path_energy_balances_foldername=snakemake.params.path_energy_balances_foldername,
        path_energy_balances_europe_filename=snakemake.params.path_energy_balances_europe_filename,
        path_energy_balances_switzerland_filename=snakemake.params.path_energy_balances_switzerland_filename,
        countries=snakemake.params.countries,
        year=snakemake.params.year,
        output_path=snakemake.output[0]
    )
