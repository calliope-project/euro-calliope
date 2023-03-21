import pandas as pd
import geopandas as gpd
import numpy as np
from shapely.geometry import Point

WGS_84 = "EPSG:4326"


def compute_emission_fractions(path_industrial_emission_data_master: str,
                               path_units: str,
                               path_mapping_unit_codes_nat_codes: str,
                               path_units_share_of_nat_demand_per_sector: str) -> None:

    """
    First, fill all NaNs and "no_matchings" in the indicator columns of df_master (such as CO2 emissions,
    NOX emissions, number of employees, ...) with 0. Then, for each tuple (unit i, sector j) compute the
    unit's share of the national pollutant emissions of that sector. This will later serve as measure for
    each units's share of the national electricity demand of a sector (demand_sectorj_uniti / demand_sectorj_national).

    Parameters:
        path_industrial_emission_data_master (str):
            File (.csv) with emission and employee data for each installation
        path_units (str):
            File (.geojson) with geometry of units
        path_mapping_unit_codes_nat_codes (str):
            File that allows mapping between unit codes and country codes
        path_units_share_of_nat_demand_per_sector (str):
            Output file with share of each unit's share of the national electricity demand of that sector
    """

    # Import units (for the respective resolution -- only necessary for all subnational resolutions!)
    # Import emissions per installation
    # Import file that allows mapping between unit codes and country codes
    units = gpd.read_file(path_units).to_crs(WGS_84).set_index("id") # or any other of the unit...geojson files.
    df_master = pd.read_csv(path_industrial_emission_data_master, index_col=0, sep=',', quotechar='"')
    mapping_unit_codes_nat_codes = pd.read_csv(path_mapping_unit_codes_nat_codes, index_col="id")

    df_master = df_master.reset_index(drop=True) #ToDo: remove.

    # Add unit code based on unit boundaries and installation coordinates
    installation_coords = gpd.GeoDataFrame(
        crs=WGS_84,
        geometry=list(map(Point, zip(df_master.pointGeometryLon, df_master.pointGeometryLat))),
        index=df_master.index)

    breakpoint()
    df_master["unitCode"] = gpd.sjoin(installation_coords, units, how="left")["index_right"]
    if (len(units.index) == 1) and (units.index[0] == "EUR"): # special case for continental level
        # In order for later aggregation over "countryCode" to work, replace countryCodes of all installations
        # within the considered unit (Note that even if the unit is the continent, not all installations necessarily
        # fall into the unit, since only selected countries can compose the continent) (i.e., installations that have
        # unitCode "EUR" not NaN)
        df_master.loc[df_master["unitCode"] == "EUR", "countryCode"] = "EUR"
    breakpoint()

    # Drop installations that do not lie in any of the units (e.g., installations in Spanish/British/French.. oversea
    # territories), installations in Iceland, installations very close to a curvy coastline (is smoothed in units
    # geometries), etc. (This affects ca. 450/68000 installations in Europe):
    df_master = df_master.dropna(subset=["unitCode"])
    # Drop installations that are misplaced in neighboring country due to curvy boarder (unit code doesnt coincide with
    # country code) (This effects ca 1 installation in Europe):
    countryCodes_from_unitCodes = (
        df_master
        .unitCode
        .apply(lambda x: mapping_unit_codes_nat_codes.loc[x, "country_code"])
    )
    df_master = df_master.drop(
        df_master[countryCodes_from_unitCodes != df_master.countryCode].index)

    # - Fill all NaNs and "no_matchings" in the indicator columns of df_master with 0 (ASSUME that no
    # emissions/employees for installations that dont report).
    # - Note that Number of Employees is an optional reporting item in the IED dataset. Can still be used as size
    # indicator if, within a nation and a sector, in each unit the same percentage of companies reports employees.
    # Otherwise apply data imputation.
    # - Group first by sector and then, within sector, group installations by country code and sum non-zero values to
    # compute national emission value of the sector
    # - Group first by sector and then, within sector, group installations by unit code and sum non-zero values to
    # compute unit's emission value of the sector
    indicatorlist = ["verified_ETS_sum", "allocatedFree_ETS_sum", "numberOfEmployees", "totalWasteQuantityTNE",
                     "ANTHRACENE", "ASANDCOMPOUNDS", "ATRAZINE", "BENZENE", "BENZO(G,H,I)PERYLENE", "CDANDCOMPOUNDS",
                     "CFCS", "CH4", "CHLORIDES", "CHLORINEANDINORGANICCOMPOUNDS", "CHLORO-ALKANES(C10-13)",
                     "CHLORPYRIFOS", "CO", "CO2", "CO2EXCLBIOMASS", "CRANDCOMPOUNDS", "CUANDCOMPOUNDS", "CYANIDES",
                     "DCE-1,2", "DCM", "DEHP", "DIURON", "ETHYLBENZENE", "ETHYLENEOXIDE", "FLUORANTHENE", "FLUORIDES",
                     "FLUORINEANDINORGANICCOMPOUNDS", "HALOGENATEDORGANICCOMPOUNDS", "HALONS", "HCB", "HCBD", "HCFCS",
                     "HCH", "HCN", "HFCS", "HGANDCOMPOUNDS", "ISOPROTURON", "LINDANE", "N2O", "NAPHTHALENE", "NH3",
                     "NIANDCOMPOUNDS", "NMVOC", "NOX", "NPANDNPES", "OCTYLPHENOLSANDOCTYLPHENOLETHOXYLATES",
                     "ORGANOTINCOMPOUNDS", "PAHS", "PBANDCOMPOUNDS", "PCBS", "PCDD+PCDF(DIOXINS+FURANS)", "PCP",
                     "PENTACHLOROBENZENE", "PER", "PFCS", "PHENOLS", "PM10", "PM2.5", "SF6", "SIMAZINE", "SOX", "TCB",
                     "TCE-1,1,1", "TCM", "TETRACHLOROETHANE-1,1,2,2", "TOC", "TOLUENE", "TOTALNITROGEN",
                     "TOTALPHOSPHORUS", "TRI", "TRIBUTYLTINANDCOMPOUNDS", "TRICHLOROMETHANE",
                     "TRIPHENYLTINANDCOMPOUNDS", "VINYLCHLORIDE", "XYLENES", "ZNANDCOMPOUNDS"]
    # ASSUME that no emissions/employees for installations that dont report:
    dfnan = pd.DataFrame(np.zeros((len(df_master.index), len(indicatorlist))), columns=indicatorlist)
    df_master = df_master.fillna(dfnan)
    # ASSUME that no emissions/employees for installations that dont report:
    functionlist = [lambda x: x.replace("no_matching", 0).astype(float).sum()] * len(indicatorlist)
    national_emissions_of_sector = (
        df_master
        .groupby(["activity_id_EUROSTAT", "countryCode"])
        .agg(dict(zip(indicatorlist, functionlist))))
    unit_emissions_of_sector = (
        df_master
        .groupby(["activity_id_EUROSTAT", "unitCode"])
        .agg(dict(zip(indicatorlist, functionlist))))
    # For each unit and sector, divide unit's emission by national emission
    fraction_unit_of_nat = unit_emissions_of_sector # to have same format
    for index, row in unit_emissions_of_sector.iterrows(): # index[0] holds activity_id_EEA, index[1] holds unitCode
        country_code = mapping_unit_codes_nat_codes.loc[index[1], "country_code"] # gives alpha 3 country code
        fraction_unit_of_nat.loc[index, :] = unit_emissions_of_sector.loc[index[0], index[1]].divide(
            national_emissions_of_sector.loc[index[0], country_code])

    # Fill nan that occur through dividing by zero (no national emissions of this pollutant)
    fraction_unit_of_nat = fraction_unit_of_nat.fillna(0)

    breakpoint()
    fraction_unit_of_nat_agg = aggregate_fractions(fraction_unit_of_nat)

    # Include missing sectors and units and assign 0 to as industrial emission fractions
    ind_sector_list = ["FC_IND_CON_E", "FC_IND_CPC_E", "FC_IND_FBT_E", "FC_IND_IS_E", "FC_IND_MAC_E", "FC_IND_MQ_E",
                       "FC_IND_NFM_E", "FC_IND_NMM_E", "FC_IND_NSP_E", "FC_IND_PPP_E", "FC_IND_TE_E", "FC_IND_TL_E",
                       "FC_IND_WP_E", "FC_OTH_AF_E", "FC_OTH_CP_E", "FC_OTH_FISH_E", "FC_OTH_HH_E", "FC_OTH_NSP_E",
                       "FC_TRA_DAVI_E", "FC_TRA_DNAVI_E", "FC_TRA_PIPE_E", "FC_TRA_RAIL_E", "FC_TRA_ROAD_E",
                       "FC_TRA_NSP_E"]
    fraction_unit_of_nat_agg = (
        fraction_unit_of_nat_agg
        + pd.DataFrame(np.zeros((len(ind_sector_list), len(units.index.tolist()))),
                       index=ind_sector_list,
                       columns=units.index.tolist()))
    fraction_unit_of_nat_agg = fraction_unit_of_nat_agg.fillna(0)

    fraction_unit_of_nat_agg.to_csv(path_units_share_of_nat_demand_per_sector)


def aggregate_fractions(fraction_unit_of_nat: pd.DataFrame) -> pd.DataFrame:

    """
    Aggregate fractions for the different pollutants into a single fraction per unit per sector.

    Parameters:
        Dataframe with unit's share of the national pollutant emissions of that sector for all pollutants.
    Returns:
        Dataframe with unit's share of the national pollutant emissions of that sector, aggregated over all pollutants.
    """

    fraction_unit_of_nat_agg = pd.DataFrame()
    for index, row in fraction_unit_of_nat.iterrows(): # index[0] holds activity_id_EEA, index[1] holds unitCode
        # Aggregation method: Only use CO2 emissions reported to ETS
        fraction_unit_of_nat_agg.loc[index] = row["verified_ETS_sum"]

    # Fill nan that occur because there is not a single installation of sector j in unit i
    fraction_unit_of_nat_agg = fraction_unit_of_nat_agg.fillna(0)

    return fraction_unit_of_nat_agg


if __name__ == "__main__":
    compute_emission_fractions(
        path_industrial_emission_data_master=snakemake.input.industrial_emission_data_master,
        path_units=snakemake.input.units,
        path_mapping_unit_codes_nat_codes=snakemake.input.mapping_unit_codes_nat_codes,
        path_units_share_of_nat_demand_per_sector=snakemake.output[0]
    )
