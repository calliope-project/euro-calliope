import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import pycountry
import re


def preprocess_emission_data(path_ied_etsidentifiers: str,
                             path_ied_general_facility_infos: str,
                             path_ied_waste: str,
                             path_ied_pollutants: str,
                             path_ied_waste_water: str,
                             path_ied_match_facil_insp_inst_insp: str,
                             reportingYear_ETS: int,
                             path_ets_installation: str,
                             path_ets_compliance: str,
                             path_mapping_nace_ets: str,
                             path_mapping_nace_eurostat: str,
                             path_industrial_emission_data_master: str,
                             path_industrial_emission_data_installation_details: str) -> None:
    """
    Preprocess IED data. Preprocess ETS data. Identify plants that occur in both datasets (matching). Join both
    datasets. Convert plants' activity codes to EUROSTAT format, convert country codes, data cleansing, reindexing.

    Parameters:
        path_ied_etsidentifiers (str):
            File needed to create matching between IED and ETS installations
        path_ied_general_facility_infos (str):
            Meta-information for IED facilities: Facility_INSPIRE_ID, pointGeometryLat, pointGeometryLon,
            NUTSRegionSourceCode, NACEMainEconomicActivityCode, mainActivityCode (for Facilities, as defined in EU's
            E-PRTR regulation), numberOfOperatingHours, numberOfEmployees
        path_ied_waste (str):
            Waste data of IED facilities: Facility_INSPIRE_ID, totalWasteQuantityTNE (reporting year = 2020). Later
            sum over wasteClassification (Hazardous, Non-hazardous) and over wasteTreatment (Disposal, Recovery).
        path_ied_pollutants (str):
            Pollutant data of IED facilities: Facility_INSPIRE_ID, pollutantCode, medium,
            totalPollutantQuantityKG (reporting year = 2020).
        path_ied_waste_water (str):
            Waste water data of IED facilities: Facility_INSPIRE_ID, pollutantCode,
            totalPollutantQuantityKg (reporting year = 2020).
        path_ied_match_facil_insp_inst_insp (str):
            File needed to relate IED facilitiy (parent) and IED installations (child): Facility_INSPIRE_ID,
            Installation_INSPIRE_ID, mainActivityCode (for Installation, as defined in EUs IED regulation). For
            matching of Facility_INSPIRE_ID and Installation_INSPIRE_ID.
        reportingYear_ETS (int):
            Year for which ETS data is used. ReportingYear for IED data is set in the rule
            manipulate_ied_etsidentifiers.
        path_ets_installation (str):
            Meta-information for ETS installations: Installation ID, Permit ID, E-PRTR ID, Lat/Lon, address, NACE
            activity code, i.a.
        path_ets_compliance (str):
            Emission data for ETS installations: verified emissions, free allocated emissions
        path_mapping_nace_ets (str):
            File needed for conversion of ETS activity codes to NACE activity codes.
        path_mapping_nace_eurostat (str):
            File needed for conversion of NACE activity codes to EUROSTAT sector codes.
        path_industrial_emission_data_master (str):
            Output file that lists all facilities reported under ETS or IED with their industrial activity code,
            Lat/Lon, ETS emissions, IED emissions, number of operating hours, number of employees, i.a.
        path_industrial_emission_data_installation_details (str):
            Supplementary output file that lists details for all ETS installations that have been matched to an IED
            facility. The background is that some IED facilities contain multiple IED installations which in turn can
            contain multiple ETS installations. Here the emissions of these installations are summed so that only a
            single value per facility is reported in the master dataset. For the individual emissions and individual
            economic activities (which can deviate from the facilities main activity) not to be lost they are reported
            in this file.
    """

    df_query1, df_master, df_query3, df_query45, df_query6 = preprocess_ied_data(path_ied_etsidentifiers,
                                                                                 path_ied_general_facility_infos,
                                                                                 path_ied_waste,
                                                                                 path_ied_pollutants,
                                                                                 path_ied_waste_water,
                                                                                 path_ied_match_facil_insp_inst_insp)
    print("Done with preprocessing IED data (Step: 1/5). For next step expect: 15 min.")
    df_ets = preprocess_ets_data(reportingYear_ETS,
                                 path_ets_installation,
                                 path_ets_compliance)
    print("Done with preprocessing ETS data (Step: 2/5). For next step expect: 5-10 min.")
    df_matching = identify_matching_ets_ied(df_query1,
                                            df_ets)
    print("Done with matching ETS and IED installations (Step: 3/5). For next step expect: 5-10 min.")
    df_master, df_master_installation_detail = join_ets_ied(df_master,
                                                            df_query3,
                                                            df_query45,
                                                            df_query6,
                                                            df_ets,
                                                            df_matching)
    print("Done with joining the datasets based on the matching (Step: 4/5). For next step expect: few sec.")
    df_master = convert_activity_codes(df_master,
                                       path_mapping_nace_ets,
                                       path_mapping_nace_eurostat)
    print("Done with conversion to EUROSTAT activity codes (Step: 5/5).")

    df_master.to_csv(path_industrial_emission_data_master, sep=',', quotechar='"')
    df_master_installation_detail.to_csv(path_industrial_emission_data_installation_details, sep=',', quotechar='"')


def preprocess_ied_data(path_ied_etsidentifiers: str,
                        path_ied_general_facility_infos: str,
                        path_ied_waste: str,
                        path_ied_pollutants: str,
                        path_ied_waste_water: str,
                        path_ied_match_facil_insp_inst_insp: str
                        ) -> "tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]":
    """
    Preprocess IED data.

    Parameters:
        path_ied_etsidentifiers (str):
            File later needed to create matching between IED and ETS installations
        path_ied_general_facility_infos (str):
            Meta-information for IED facilities: Facility_INSPIRE_ID, pointGeometryLat, pointGeometryLon,
            NUTSRegionSourceCode, NACEMainEconomicActivityCode, mainActivityCode (for Facilities, as defined in EU's
            E-PRTR regulation), numberOfOperatingHours, numberOfEmployees
        path_ied_waste (str):
            Waste data of IED facilities: Facility_INSPIRE_ID, totalWasteQuantityTNE (reporting year = 2020). Later
            sum over wasteClassification (Hazardous, Non-hazardous) and over wasteTreatment (Disposal, Recovery).
        path_ied_pollutants (str):
            Pollutant data of IED facilities: Facility_INSPIRE_ID, pollutantCode, medium,
            totalPollutantQuantityKG (reporting year = 2020).
        path_ied_waste_water (str):
            Waste water data of IED facilities: Facility_INSPIRE_ID, pollutantCode,
            totalPollutantQuantityKg (reporting year = 2020).
        path_ied_match_facil_insp_inst_insp (str):
            File needed to relate IED facilitiy (parent) and IED installations (child): Facility_INSPIRE_ID,
            Installation_INSPIRE_ID, mainActivityCode (for Installation, as defined in EUs IED regulation). For
            matching of Facility_INSPIRE_ID and Installation_INSPIRE_ID.

    Returns:
        df_query1 (pd.DataFrame):
            DataFrame with preprocessed ETS Identifiers as reported in the IED dataset to link IED installations to
            ETS installations.
        df_master (pd.DataFrame):
            DataFrame with preprocessed meta-information of IED facilities. Information from other DataFrames are
            added to this one.
        df_query3 (pd.DataFrame):
            DataFrame with preprocessed waste data for IED facilities.
        df_query45 (pd.DataFrame):
            DataFrame with preprocessed pollutant and waste water data for IED facilities.
        df_query6 (pd.DataFrame):
            DataFrame with preprocessed information to relate IED facilitiy (parent) and IED installations (child).
    """

    # Further processing QUERY 1
    # Drop duplicate rows: drop if all fields (besides internal numeration identifier) are equal
    df_query1 = pd.read_csv(path_ied_etsidentifiers, sep=';', quotechar='"')
    df_query1 = df_query1.drop_duplicates(subset=['Installation_INSPIRE_ID', 'reportingYear', 'ETSIdentifier',
                                                  'pointGeometryLat', 'pointGeometryLon', 'nameOfFeature',
                                                  'streetName', 'buildingNumber', 'city', 'postalCode',
                                                  'countryCode', 'conv_ETSIdentifier'])

    # Further processing QUERY 2
    # The dataframe is named master since the infos from all other tables will gradually be added to this table.
    df_master = pd.read_csv(path_ied_general_facility_infos, index_col=["Facility_INSPIRE_ID"], sep=';', quotechar='"')
    # Rename mainActivityCode to mainActivityCode_IED_Faci to avoid confusion with mainActivityCode_IED_Inst and
    # activity_id_ETS and nace_id_ETS and NACEMainEconomicActivityCode
    df_master.rename(columns={'mainActivityCode': 'mainActivityCode_IED_Faci'}, inplace=True)
    # A few facilities (<20 facilities) occur multiple times in this table since they have multiple
    # NACEMainActivityCodes. However it is hardly possible to determine which is the main activity of a facility,
    # since e.g. the emissions are reported for the entire facility, not per activity. Therefore, we group the
    # facilities by Facility_INSPIRE_ID and keep the first NACE activity code. Alternatively, the "mainActivityCode"
    # (only one of these is provided per facility) could be used to identify the correct NACEMainActivityCode.
    df_master = (
        df_master
        .groupby(['Facility_INSPIRE_ID'], as_index=True)
        .agg({"NACEMainEconomicActivityCode": 'first', "mainActivityCode_IED_Faci": 'first',
              "pointGeometryLat": 'first', "pointGeometryLon": 'first', "NUTSRegionSourceCode": 'first',
              "countryCode": 'first', "reportingYear": 'first', "numberOfOperatingHours": 'first',
              "numberOfEmployees": 'first'})
    )

    # Further processing QUERY 3
    df_query3 = pd.read_csv(path_ied_waste, sep=';', quotechar='"')
    # Facilities report hazardous and non-hazardous waste, as well as disposed and recycled waste. Here, all wastes
    # are summed for each facility.
    df_query3 = df_query3.groupby(['Facility_INSPIRE_ID'], as_index=False).agg({'totalWasteQuantityTNE': 'sum'})
    # Set index column
    df_query3 = df_query3.set_index("Facility_INSPIRE_ID")

    # Further processing QUERY 4
    df_query4 = pd.read_csv(path_ied_pollutants, sep=';', quotechar='"')
    df_query4 = df_query4.drop(columns="reportingYear")

    # Further processing QUERY 5
    df_query5 = pd.read_csv(path_ied_waste_water, sep=';', quotechar='"')
    # Prepare for merge with IED_Pollutants by assigning the lable WASTEWATER
    df_query5["medium"] = "WASTEWATER"
    # Drop info of reporting year, no longer needed
    df_query5 = df_query5.drop(columns="reportingYear")

    # Further processing QUERY 6
    # Load and set index column
    df_query6 = (
        pd.read_csv(path_ied_match_facil_insp_inst_insp, index_col=["Installation_INSPIRE_ID"], sep=';', quotechar='"')
    )

    # Merge data from QUERY 4 (Pollutants released to Air, Land, Water) and QUERY 5 (Wastewater)
    df_query45 = pd.concat([df_query5, df_query4])
    # For each pollutant sum over emissions to Air, Land, Water, Wastewater. Set index column.
    df_query45 = df_query45.set_index(["Facility_INSPIRE_ID", "pollutantCode", "medium"])
    df_query45 = df_query45.groupby(['Facility_INSPIRE_ID', 'pollutantCode']).agg({'totalPollutantQuantityKg': 'sum'})
    df_query45 = df_query45["totalPollutantQuantityKg"].unstack(level=-1)

    return df_query1, df_master, df_query3, df_query45, df_query6


def preprocess_ets_data(reportingYear_ETS: int, path_ets_installation: str, path_ets_compliance: str) -> pd.DataFrame:
    """
    Preprocess ETS data.

    First step: The ETS dataset (file "installation.csv") provides ETS installation IDs in the format DE_245. Extend
        these to ISO code (DE000000000000245).
    Second step: Add the emission data (Fields "verified" (verified emissions) and "allocatedFree" (freely allocated
        emissions) in "compliance.csv") to the installations in "installation.csv" for a given year. The two tables
        can be matched as follows: field "id" in "installation.csv" corrsponds to field "installation_id" in
        "compliance.csv". The fields "verified" and "allocatedFree" remain empty for installations that did not report
        any values in the chosen year.
    Note: The current implementation of the second step is rather slow (expect 15 min).

    Parameters:
        reportingYear_ETS (int):
            Year for which ETS data is used. ReportingYear for IED data is set in the rule
            manipulate_ied_etsidentifiers.
        path_ets_installation (str):
            Meta-information for ETS installations: Installation ID, Permit ID, E-PRTR ID, Lat/Lon, address, NACE
            activity code, i.a.
        path_ets_compliance (str):
            Emission data for ETS installations: verified emissions, free allocated emissions

    Returns:
        df_inst (pd.DataFrame):
            ETS dataset as DataFrame. Includes the fields ISO ID, Permit ID, and E-PRTR ID - the three ID types that
            might be reported as ETS Identifiers in the IED dataset.
    """

    df_inst = pd.read_csv(path_ets_installation, sep=',', quotechar='"')
    df_compl = pd.read_csv(path_ets_compliance, sep=',', quotechar='"')

    df_inst.insert(len(df_inst.columns), "id_ISO", "")
    df_inst.insert(len(df_inst.columns), "verified", "")
    df_inst.insert(len(df_inst.columns), "allocatedFree", "")

    for index, row in df_inst.iterrows():
        # Add column with ISO installation ID
        numeral = re.sub("\D", "", row["id"])
        ISO_code = row["registry_id"] + ("0" * (15 - len(numeral))) + numeral
        df_inst.loc[index, "id_ISO"] = ISO_code

        # Add columns with emission data
        compl_rowindex = (
            df_compl
            .index[(df_compl["installation_id"] == row["id"])
                   & (df_compl["year"] == reportingYear_ETS)
                   & (df_compl["reportedInSystem"] == "EUETS")]
            .tolist()
        )
        df_inst.loc[index, "verified"] = df_compl.loc[compl_rowindex[0], "verified"]
        df_inst.loc[index, "allocatedFree"] = df_compl.loc[compl_rowindex[0], "allocatedFree"]

    return df_inst


def identify_matching_ets_ied(df_query1: pd.DataFrame, df_ets: pd.DataFrame) -> pd.DataFrame:
    """
    Identify plants that occur in both datasets (matching).

    First identify potential matchings via comparison of conv_ETSIdentifier (IED dataset) and id_ISO, permitID, and
        eprtrID (ETS dataset)
    Then verify matching by comparison of postal codes and geographic coordinates.
    In the resulting file, the Installation_INSPIRE_ID is the unique identifier in the IED dataset and the id_ISO is
        the unique identifier in the ETS dataset.
    Note: Current implementation takes ca. 6 minutes for execution.

    Parameters:
        df_query1 (pd.DataFrame):
            DataFrame with preprocessed ETS Identifiers as reported in the IED dataset to link IED installations to
            ETS installations.
        df_ets (pd.DataFrame):
            Preprocessed ETS dataset as DataFrame. Includes the fields ISO ID, Permit ID, and E-PRTR ID - the three
            ID types that might be reported as ETS Identifiers in the IED dataset.

    Returns:
        df_matching (pd.DataFrame):
            DataFrame that relates ETS installations (via their ISO ID) to IED facilities (via their
            Facility_INSPIRE_ID).
    """

    df_ied = df_query1

    # Create dataframe to record matches
    df_matching = pd.DataFrame()

    for index, row in df_ied.iterrows():
        conv_ETSIdentifier = row["conv_ETSIdentifier"]

        # Skip installations with conv_ETSIdentifier = 0, since 0 is not a valid Identifier in ETS dataset anyway.
        # Also the list_of_matching_indexes would be very long for these installations.
        if conv_ETSIdentifier == '0':
            continue

        list_of_matching_indexes = (
            df_ets
            .index[((df_ets["id_ISO"] == conv_ETSIdentifier)
                   | (df_ets["permitID"] == conv_ETSIdentifier)
                   | (df_ets["eprtrID"] == conv_ETSIdentifier))]
            .tolist()
        )

        for ets_index in list_of_matching_indexes:

            # Set default values
            matched = False
            matchedISO = ""
            verifLatLong, verifPostalCode = [False, False]
            viaISO, viaPermitID, viaEprtrID = [False, False, False]

            # Sometimes multiple IED installations map to the same ETS account and therefore emission value. Thus,
            # only the first matching of each ISO_ID is allowed, to avoid double consideration. This is a time
            # consuming step.
            if not df_matching.empty:
                if df_ets.loc[ets_index, "id_ISO"] in df_matching["id_ISO"].values:
                    continue

            # Verify matches based on coordinates:
            geometry_ets = [Point(df_ets.loc[ets_index, "longitudeGoogle"], df_ets.loc[ets_index, "latitudeGoogle"])]
            points_df_ets = gpd.GeoDataFrame({'geometry': geometry_ets}, crs='EPSG:4326').to_crs('EPSG:3035')
            geometry_ied = [Point(row["pointGeometryLon"], row["pointGeometryLat"])]
            points_df_ied = gpd.GeoDataFrame({'geometry': geometry_ied}, crs='EPSG:4326').to_crs('EPSG:3035')
            distance = points_df_ets.distance(points_df_ied)
            if distance[0] < 2000:
                verifLatLong = True

            # Verify matches based on postalCode:
            if df_ied.loc[index, "postalCode"] == df_ets.loc[ets_index, "postalCode"]:
                verifPostalCode = True

            # If matching verfied, replace default values with matching infos:
            if verifLatLong or verifPostalCode:
                matched = True
                matchedISO = df_ets.loc[ets_index, "id_ISO"]

                # Report via which ID the installations were matched.
                if row["conv_ETSIdentifier"] == df_ets.loc[ets_index, "id_ISO"]:
                    viaISO = True
                if row["conv_ETSIdentifier"] == df_ets.loc[ets_index, "permitID"]:
                    viaPermitID = True
                if row["conv_ETSIdentifier"] == df_ets.loc[ets_index, "eprtrID"]:
                    viaEprtrID = True

                list_with_match = [row["Installation_INSPIRE_ID"], row["conv_ETSIdentifier"], matched, matchedISO,
                                   viaISO, viaPermitID, viaEprtrID, verifLatLong, verifPostalCode]
                df_matching = (
                    df_matching
                    .append(pd.DataFrame([list_with_match], columns=["Installation_INSPIRE_ID", "conv_ETSIdentifier",
                                                                     "matched", "id_ISO", "viaISO", "viaPermitID",
                                                                     "viaEprtrID", "verifLatLong", "verifPostalCode"]))
                )

    return df_matching


def join_ets_ied(df_master: pd.DataFrame, df_query3: pd.DataFrame, df_query45: pd.DataFrame, df_query6: pd.DataFrame,
                 df_ets: pd.DataFrame, df_matching: pd.DataFrame) -> "tuple[pd.DataFrame, pd.DataFrame]":
    """
    Join ETS and IED dataset.

    Script creates master dataframe (df_master) that lists size indicators as well as meta information for all
        reporting plants from both datasets.
    Size indicators of a plant include: verified CO2 emissions (ETS dataset), freely allocated certificates (ETS
        dataset), number of employees (IED), total waste in tons (IED), emissions of 76 pollutants in kg (IED).
    Meta information of a plant include: Facility_INSPIRE_ID (unique IED ID), id_ISO (unique ETS ID),
        NACEMainEconomicActivityCode, mainActivityCode_IED_Inst, mainActivityCode_IED_Faci, activity_id_ETS,
        pointGeometryLat, pointGeometryLon, NUTSRegionSourceCode, countryCode, numberOfOperatingHours
    For facilities from the IED dataset, that could not be matched with an ETS installation the following fields
        remain empty (i.e., "no_matching"): id_ISO, mainActivityCode_IED_Inst, activity_id_ETS, verified_ETS,
        allocatedFree_ETS
    For installations from the ETS dataset, that could not be matched with an IED facility the following fields remain
        empty (i.e., "no_matching"): Facility_INSPIRE_ID, mainActivityCode_IED_Inst, mainActivityCode_IED_Faci,
        NUTSRegionSourceCode, numberOfOperatingHours, numberOfEmployees, totalWasteQuantity, pollutant emissions
    Note that "no_matching" marks those fields that remain empty because no matching installation from the other
        dataset could be identified and therefore the information for this facility is incomplete. On the contrary,
        other missing values, NaNs, and 0 stem from the original dataset and can represent either real zeros or
        unreported data points.
    Note that CO2 emissions are reported in both, the IED and the ETS dataset. However, typically plants that must
        report their CO2 emissions under ETS regulations must not report them under IED regulations and vice versa. If
        a facility reports both, the emissions stem from different subparts of the facility: e.g., one obliged to
        report to ETS and one obliged to report to IED. No risk of double consideration.
    Note that multiple IED installations may exist per IED facility, and that multiple ETS reports may exist for the
        same IED installation. In these cases the field entries are arranged in lists. emissions are summed over all
        installations for each facility.
    An additional dataframe (df_master_installation_detail) is created that records and preserves the information of
        all installations (e.g., IED installations appear twice if they have two ETS records assigned), also providing
        their parent facility ID.
    Note: The current implementation of this script takes 5-10 minutes.

    Parameters:
        df_master (pd.DataFrame):
            DataFrame with preprocessed meta-information of IED facilities. Information from other DataFrames are
            added to this one.
        df_query3 (pd.DataFrame):
            DataFrame with preprocessed waste data for IED facilities.
        df_query45 (pd.DataFrame):
            DataFrame with preprocessed pollutant and waste water data for IED facilities.
        df_query6 (pd.DataFrame):
            DataFrame with preprocessed information to relate IED facilitiy (parent) and IED installations (child).
        df_ets (pd.DataFrame):
            DataFrame with preprocessed meta and emission information on installations reporting under the ETS
            regulation.
        df_matching (pd.DataFrame):
            DataFrame that relates ETS installations (via their ISO ID) to IED facilities (via their
            Facility_INSPIRE_ID).

    Returns:
        df_master (pd.DataFrame):
            DataFrame that lists all facilities reported under ETS or IED with their industrial activity code,
            Lat/Lon, ETS emissions, IED emissions, number of operating hours, number of employees, i.a.
        df_master_installation_detail (pd.DataFrame):
            DataFrame that lists details for all ETS installations that have been matched to an IED facility. The
            background is that some IED facilities contain multiple IED installations which in turn can contain
            multiple ETS installations. Here the emissions of these installations are summed so that only a single
            value per facility is reported in the master dataset. For the individual emissions and individual economic
            activities (which can deviate from the facilities main activity) not to be lost they are reported in this
            file.
    """

    # Set index in ETS data.
    df_ets = df_ets.set_index("id_ISO")

    # Merge general installation infos with amount of waste and amount of pollutants
    df_master = df_master.merge(df_query3, how='outer', on='Facility_INSPIRE_ID')
    df_master = df_master.merge(df_query45, how='outer', on='Facility_INSPIRE_ID')

    # Initialize df_master_installation_detail
    df_master_installation_detail = pd.DataFrame()

    # df_master contains all IED facilities in the reporting year, now information from the ETS dataset is to be added.
    # First create columns in df_master to hold information specific to the ETS dataset (e.g., verified emissions).
    # Set "no_matching" as a default value which will be overriden for IED facilities that match with an ETS facility.
    # Then iterate through matched IED installations, and group them with their parent IED facility.

    df_master.insert(0, "id_ISO", "no_matching")
    df_master.insert(2, "mainActivityCode_IED_Inst", "no_matching")
    df_master.insert(4, "activity_id_ETS", "no_matching")
    df_master.insert(9, "verified_ETS", "no_matching")
    df_master.insert(10, "allocatedFree_ETS", "no_matching")
    df_master.drop("reportingYear", axis=1, inplace=True)

    for index, row in df_matching.iterrows():

        # Gather information about matched installation for df_master and df_master_installation_detail
        facility_INSPIRE_ID = df_query6.loc[row["Installation_INSPIRE_ID"], "Parent_Facility_INSPIRE_ID"]
        mainActivityCode_IED_Inst = df_query6.loc[row["Installation_INSPIRE_ID"], "mainActivityCode"]
        id_ISO = row["id_ISO"]
        activity_id_ETS = df_ets.loc[id_ISO, "activity_id"]
        verified_ETS = df_ets.loc[id_ISO, "verified"]
        allocatedFree_ETS = df_ets.loc[id_ISO, "allocatedFree"]

        # Gather information about matched installation for df_master_installation_detail exclusively
        nace_id_ETS = df_ets.loc[id_ISO, "nace_id"]
        latitudeGoogle_ETS = df_ets.loc[id_ISO, "latitudeGoogle"]
        longitudeGoogle_ETS = df_ets.loc[id_ISO, "longitudeGoogle"]
        NACEMainEconomicActivityCode_IED = df_master.loc[facility_INSPIRE_ID, "NACEMainEconomicActivityCode"]
        pointGeometryLat_IED = df_master.loc[facility_INSPIRE_ID, "pointGeometryLat"]
        pointGeometryLon_IED = df_master.loc[facility_INSPIRE_ID, "pointGeometryLon"]

        # Write to df_master. If the first installation is added to a facility, replace "no_matching" with empty list,
        # that is appended in the following.
        if df_master.loc[facility_INSPIRE_ID, "id_ISO"] == str("no_matching"):
            df_master.at[facility_INSPIRE_ID, "id_ISO"] = list()
            df_master.at[facility_INSPIRE_ID, "mainActivityCode_IED_Inst"] = list()
            df_master.at[facility_INSPIRE_ID, "activity_id_ETS"] = list()
            df_master.at[facility_INSPIRE_ID, "verified_ETS"] = list()
            df_master.at[facility_INSPIRE_ID, "allocatedFree_ETS"] = list()
        df_master.at[facility_INSPIRE_ID, "id_ISO"] = df_master.loc[facility_INSPIRE_ID, "id_ISO"] + [id_ISO]
        df_master.at[facility_INSPIRE_ID, "mainActivityCode_IED_Inst"] = (
            df_master.loc[facility_INSPIRE_ID, "mainActivityCode_IED_Inst"] + [mainActivityCode_IED_Inst]
        )
        df_master.at[facility_INSPIRE_ID, "activity_id_ETS"] = (
            df_master.loc[facility_INSPIRE_ID, "activity_id_ETS"] + [activity_id_ETS]
        )
        df_master.at[facility_INSPIRE_ID, "verified_ETS"] = (
            df_master.loc[facility_INSPIRE_ID, "verified_ETS"] + [verified_ETS]
        )
        df_master.at[facility_INSPIRE_ID, "allocatedFree_ETS"] = (
            df_master.loc[facility_INSPIRE_ID, "allocatedFree_ETS"] + [allocatedFree_ETS]
        )

        # Write to df_master_installation_detail
        df_master_installation_detail = (
            pd
            .concat([df_master_installation_detail,
                     pd.DataFrame([[facility_INSPIRE_ID, id_ISO, pointGeometryLat_IED, pointGeometryLon_IED,
                                    latitudeGoogle_ETS, longitudeGoogle_ETS, NACEMainEconomicActivityCode_IED,
                                    nace_id_ETS, mainActivityCode_IED_Inst, activity_id_ETS, verified_ETS,
                                    allocatedFree_ETS]
                                   ], columns=["Facility_INSPIRE_ID", "id_ISO", "pointGeometryLat_IED",
                                               "pointGeometryLon_IED", "latitudeGoogle_ETS", "longitudeGoogle_ETS",
                                               "NACEMainEconomicActivityCode_IED", "nace_id_ETS",
                                               "mainActivityCode_IED", "activity_id_ETS", "verified_ETS",
                                               "allocatedFree_ETS"]
                                  )
                     ])
        )

    # Compute sum for verified emissions and freely allocated certificates in df_master (if one IED facility has
    # multiple IED or ETS installations)
    df_master.insert(
        10,
        "verified_ETS_sum",
        df_master["verified_ETS"].apply(lambda x: sum(x) if x != str("no_matching") else "no_matching")
    )
    df_master.insert(
        12,
        "allocatedFree_ETS_sum",
        df_master["allocatedFree_ETS"].apply(lambda x: sum(x) if x != str("no_matching") else "no_matching")
    )

    # Set correct indexes
    df_master.insert(0, "Facility_INSPIRE_ID", df_master.index) # Make facility_INSPIRE_ID a non-index column again
    df_master_installation_detail = df_master_installation_detail.reset_index(drop=True)

    # Add ETS installations that dont have match with facilities
    for index, row in df_ets.iterrows():

        if index in df_matching["id_ISO"].values:
            continue
        facility_INSPIRE_ID = "no_matching"
        id_ISO = index
        nace_id_ETS = row["nace_id"]
        latitudeGoogle_ETS = row["latitudeGoogle"]
        longitudeGoogle_ETS = row["longitudeGoogle"]
        countryCode = row["country_id"]
        activity_id_ETS = df_ets.loc[id_ISO, "activity_id"]
        mainActivityCode_IED_Faci = "no_matching"
        mainActivityCode_IED_Inst = "no_matching"
        NUTSRegionSourceCode = "no_matching"
        verified_ETS = df_ets.loc[id_ISO, "verified"]
        allocatedFree_ETS = df_ets.loc[id_ISO, "allocatedFree"]

        columns = ["Facility_INSPIRE_ID", "id_ISO", "pointGeometryLat", "pointGeometryLon",
                   "NACEMainEconomicActivityCode", "mainActivityCode_IED_Faci", "mainActivityCode_IED_Inst",
                   "activity_id_ETS", "NUTSRegionSourceCode", "countryCode", "verified_ETS", "allocatedFree_ETS",
                   "verified_ETS_sum", "allocatedFree_ETS_sum"]
        df_master = (
            pd
            .concat([df_master,
                     pd.DataFrame([[facility_INSPIRE_ID, id_ISO, latitudeGoogle_ETS, longitudeGoogle_ETS, nace_id_ETS,
                                    mainActivityCode_IED_Faci, mainActivityCode_IED_Inst, activity_id_ETS,
                                    NUTSRegionSourceCode, countryCode, [verified_ETS], [allocatedFree_ETS],
                                    verified_ETS, allocatedFree_ETS]
                                   ], columns=columns)])
        )

    # Drop Facility_Inspire_ID as index (since only an index to the IED facilities) and create integer index.
    df_master = df_master.reset_index(drop=True)

    # For ETS installations without matching with IED facility, add "no_matching" for fields that come from IED dataset
    df_master.loc[df_master.index[(df_master["Facility_INSPIRE_ID"] == "no_matching")].tolist(),
                  "numberOfOperatingHours":"ZNANDCOMPOUNDS"] = "no_matching"

    return df_master, df_master_installation_detail


def convert_activity_codes(df_master: pd.DataFrame,
                           path_mapping_nace_ets: str,
                           path_mapping_nace_eurostat: str) -> pd.DataFrame:
    """
    Convert plants' activity codes to EUROSTAT format. Convert country code from alpha 2 to alpha 3. Remove
    installations with no geo information and reset index.

    In the Master dataset, some unmatched ETS installations do not have NACE activity codes. However, they do have
    ETS's activity_IDs. This script translates and adds the activity codes to NACE IDs using the mapping made
    available via the csv file. In a second step, for all installations the NACE codes are translated into the
    categories used by the EUROSTAT in their energy balances. Note: Create the csv in UTF-8 format so that leading
    and trailing zeros are not cut.

    Parameters:
        df_master (pd.DataFrame):
            DataFrame that lists all facilities reported under ETS or IED with their industrial activity code,
            Lat/Lon, ETS emissions, IED emissions, number of operating hours, number of employees, i.a.
        path_mapping_nace_ets (pd.DataFrame):
            Path to csv file with mapping between NACE codes and activity codes as reported in the ETS dataset.
        path_mapping_nace_eurostat (pd.DataFrame):
            Path to csv file with mapping between NACE codes and activity codes as chosen by EUROSTAT in their
            energy balances.

    Returns:
        df_master (pd.DataFrame):
            Same as input argument. But with the industrial activity code also reported in the format used by EUROSTAT.
    """

    df_mapping_NACE_ETS = pd.read_csv(path_mapping_nace_ets, sep=';', quotechar='"')
    df_mapping_NACE_EUROSTAT = pd.read_csv(path_mapping_nace_eurostat, sep=';', quotechar='"')
    # Delete leading zeroes, keep trailing zeros. Now has Master's format.
    df_mapping_NACE_EUROSTAT["nace_id"] = df_mapping_NACE_EUROSTAT["nace_id"].apply(lambda x: x.lstrip("0"))

    df_master.insert(3, "activity_id_EUROSTAT", "")

    for index, row in df_master.iterrows():
        if str(row["NACEMainEconomicActivityCode"]) == "nan":
            NACE_ETS_rowindex = (
                df_mapping_NACE_ETS
                .index[df_mapping_NACE_ETS["activity_id_ETS"] == int(row["activity_id_ETS"])]
                .tolist()
            )
            df_master.loc[index, "NACEMainEconomicActivityCode"] = (
                str(df_mapping_NACE_ETS.loc[NACE_ETS_rowindex[0], "nace_id"])
            )

    for index, row in df_master.iterrows():
        NACE_EUROSTAT_rowindex = (
            df_mapping_NACE_EUROSTAT
            .index[df_mapping_NACE_EUROSTAT["nace_id"] == str(row["NACEMainEconomicActivityCode"])]
            .tolist()
        )
        df_master.at[index, "activity_id_EUROSTAT"] = (
            df_mapping_NACE_EUROSTAT.loc[NACE_EUROSTAT_rowindex[0], "EUROSTAT_EB_category"]
        )

    # Data cleaning:
    # Remove "EU" row, convert countryCodes to alpha3, convert countryCode XI (Northern Ireland) to GBR
    df_master = df_master[df_master.countryCode != "EU"]
    df_master["countryCode"] = (
        df_master["countryCode"]
        .apply(lambda iso2: "GBR" if iso2 == "XI" else pycountry.countries.lookup(iso2).alpha_3)
    )
    # Remove installations with no geographic coordinates, as required for spatial joint in
    # compute_emissions_fractions.
    df_master = df_master.dropna(axis=0, subset=['pointGeometryLon'])
    df_master = df_master.dropna(axis=0, subset=['pointGeometryLat'])
    # Reset index
    df_master.reset_index(drop=True)

    return df_master


if __name__ == "__main__":
    preprocess_emission_data(
        path_ied_etsidentifiers=snakemake.input.ied_etsidentifiers,
        path_ied_general_facility_infos=snakemake.input.ied_general_facility_infos,
        path_ied_waste=snakemake.input.ied_waste,
        path_ied_pollutants=snakemake.input.ied_pollutants,
        path_ied_waste_water=snakemake.input.ied_waste_water,
        path_ied_match_facil_insp_inst_insp=snakemake.input.ied_match_facil_insp_inst_insp,
        reportingYear_ETS=snakemake.params.reportingYear_ETS,
        path_ets_installation=snakemake.input.ets_installation,
        path_ets_compliance=snakemake.input.ets_compliance,
        path_mapping_nace_ets=snakemake.params.mapping_nace_ets,
        path_mapping_nace_eurostat=snakemake.params.mapping_nace_eurostat,
        path_industrial_emission_data_master=snakemake.output.industrial_emission_data_master,
        path_industrial_emission_data_installation_details=(
            snakemake.output.industrial_emission_data_installation_details
        )
    )
