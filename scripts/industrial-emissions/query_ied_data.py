import pypyodbc
import csv


def query_ied_data(path_ied_database: str, path_ied_etsidentifiers: str, path_ied_general_facility_infos: str,
                   path_ied_waste: str, path_ied_pollutants: str, path_ied_waste_water: str,
                   path_ied_match_facil_insp_inst_insp: str) -> None:
    """
    For loading pollution data from the IED database (.accdb format), performing queries, and exporting them to .csv.
    The following queries are performed using pypyodbc. Instead, they could also be performed manually in MS Access,
        if access to the software is available. However, when exporting manually, in the export file all numbers will
        be truncated to two decimal digits.
    Microsoft Access and also pypyodbc's Microsoft Access Driver only work for Windows. Therefore this rule can only
        be executed on Windows machines.

    Perform the following queries:
    QUERY 1: (IED_ETSIdentifiers) -------------- Extract the table ETSIdentifiers (links IED installations to ETS
                                                 installations) and enrich it with the geographic info for the IED
                                                 installations, later used for match verification.
    QUERY 2: (IED_General_Facility_Infos) ------ Facility_INSPIRE_ID, pointGeometryLat, pointGeometryLon,
                                                 NUTSRegionSourceCode, NACEMainEconomicActivityCode, mainActivityCode
                                                 (for Facilities, as defined in EU's E-PRTR regulation),
                                                 numberOfOperatingHours, numberOfEmployees
    QUERY 3: (IED_Waste) ----------------------- Facility_INSPIRE_ID, totalWasteQuantityTNE (reporting year = 2020).
                                                 Later sum over wasteClassification (Hazardous, Non-hazardous) and
                                                 over wasteTreatment (Disposal, Recovery).
    QUERY 4: (IED_Pollutants) ------------------ Facility_INSPIRE_ID, pollutantCode, medium, totalPollutantQuantityKG
                                                 (reporting year = 2020).
    QUERY 5: (IED_Waste_Water) ----------------- Facility_INSPIRE_ID, pollutantCode, totalPollutantQuantityKg
                                                 (reporting year = 2020).
    QUERY 6: (IED_Match_Facil_INSP_Inst_INSP) -- Facility_INSPIRE_ID, Installation_INSPIRE_ID, mainActivityCode (for
                                                 Installation, as defined in EUs IED regulation). For matching of
                                                 Facility_INSPIRE_ID and Installation_INSPIRE_ID.

    Parameters:
        path_ied_database (str):
            Unzipped .accdb file with IED database
        path_ied_etsidentifiers (str):
            Output csv file of query 1
        path_ied_general_facility_infos (str):
            Output csv file of query 2
        path_ied_waste (str):
            Output csv file of query 3
        path_ied_pollutants (str):
            Output csv file of query 4
        path_ied_waste_water (str):
            Output csv file of query 5
        path_ied_match_facil_insp_inst_insp (str):
            Output csv file of query 6
    """

    # Specify the year for which the IED data shall be extracted. Note that rule "manipulate_ets_identifiers" needs to
    # be adapted if dataset from other year is used.
    reportingYear = 2020

    # Establish MS ACCESS database connection
    pypyodbc.lowercase = False
    conn = pypyodbc.connect(
        r"Driver={Microsoft Access Driver (*.mdb, *.accdb)};" + "Dbq=" + path_ied_database + ";")

    # QUERY 1:
    # Open cursor and execute SQL. Then open csv and iterate through results. Then close cursor.
    cur = conn.cursor()
    cur.execute("SELECT [3f2_ETSIdentifiers].ETSIdentifierId, [3f2_ETSIdentifiers].Installation_INSPIRE_ID, "
                "[3f2_ETSIdentifiers].reportingYear, [3f2_ETSIdentifiers].ETSIdentifier, "
                "[2_ProductionFacility].pointGeometryLat, [2_ProductionFacility].pointGeometryLon, "
                "[2_ProductionFacility].nameOfFeature, [2_ProductionFacility].streetName, "
                "[2_ProductionFacility].buildingNumber, [2_ProductionFacility].city, "
                "[2_ProductionFacility].postalCode, [2_ProductionFacility].countryCode FROM (2_ProductionFacility "
                "INNER JOIN 3_ProductionInstallation ON [2_ProductionFacility].Facility_INSPIRE_ID = "
                "[3_ProductionInstallation].Parent_Facility_INSPIRE_ID) INNER JOIN 3f2_ETSIdentifiers ON "
                "[3_ProductionInstallation].Installation_INSPIRE_ID = [3f2_ETSIdentifiers].Installation_INSPIRE_ID;")
    with open(path_ied_etsidentifiers, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # QUERY 2:
    cur = conn.cursor()
    cur.execute("SELECT [2_ProductionFacility].Facility_INSPIRE_ID, [2_ProductionFacility].pointGeometryLat, "
                "[2_ProductionFacility].pointGeometryLon, [2_ProductionFacility].NUTSRegionSourceCode, "
                "[2_ProductionFacility].countryCode, [2c_Function].NACEMainEconomicActivityCode, "
                "[2_ProductionFacility].mainActivityCode, [2a_ProductionFacilityDetails].reportingYear, "
                "[2a_ProductionFacilityDetails].numberOfOperatingHours, "
                "[2a_ProductionFacilityDetails].numberOfEmployees FROM (2_ProductionFacility INNER JOIN "
                "2a_ProductionFacilityDetails ON [2_ProductionFacility].Facility_INSPIRE_ID = "
                "[2a_ProductionFacilityDetails].Facility_INSPIRE_ID) INNER JOIN 2c_Function ON "
                "[2_ProductionFacility].Facility_INSPIRE_ID = [2c_Function].Facility_INSPIRE_ID WHERE "
                "((([2a_ProductionFacilityDetails].reportingYear)=" + str(reportingYear) + "));")
    with open(path_ied_general_facility_infos, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # QUERY 3:
    cur = conn.cursor()
    cur.execute("SELECT [2_ProductionFacility].Facility_INSPIRE_ID, [2h_OffsiteWasteTransfer].reportingYear, "
                "[2h_OffsiteWasteTransfer].totalWasteQuantityTNE FROM 2_ProductionFacility INNER JOIN "
                "2h_OffsiteWasteTransfer ON [2_ProductionFacility].Facility_INSPIRE_ID = "
                "[2h_OffsiteWasteTransfer].Facility_INSPIRE_ID WHERE "
                "((([2h_OffsiteWasteTransfer].reportingYear)=" + str(reportingYear) + "));")
    with open(path_ied_waste, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # QUERY 4:
    cur = conn.cursor()
    cur.execute("SELECT [2_ProductionFacility].Facility_INSPIRE_ID, [2f_PollutantRelease].reportingYear, "
                "[2f_PollutantRelease].pollutantCode, [2f_PollutantRelease].medium, "
                "[2f_PollutantRelease].totalPollutantQuantityKg FROM 2_ProductionFacility INNER JOIN "
                "2f_PollutantRelease ON [2_ProductionFacility].Facility_INSPIRE_ID = "
                "[2f_PollutantRelease].Facility_INSPIRE_ID WHERE "
                "((([2f_PollutantRelease].reportingYear)=" + str(reportingYear) + "));")
    with open(path_ied_pollutants, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # QUERY 5:
    cur = conn.cursor()
    cur.execute("SELECT [2_ProductionFacility].Facility_INSPIRE_ID, [2g_OffsitePollutantTransfer].reportingYear, "
                "[2g_OffsitePollutantTransfer].pollutantCode, [2g_OffsitePollutantTransfer].totalPollutantQuantityKg "
                "FROM 2_ProductionFacility INNER JOIN 2g_OffsitePollutantTransfer ON "
                "[2_ProductionFacility].Facility_INSPIRE_ID = [2g_OffsitePollutantTransfer].Facility_INSPIRE_ID WHERE "
                "((([2g_OffsitePollutantTransfer].reportingYear)=" + str(reportingYear) + "));")
    with open(path_ied_waste_water, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # QUERY 6:
    cur = conn.cursor()
    cur.execute("SELECT [3_ProductionInstallation].Parent_Facility_INSPIRE_ID, "
                "[3_ProductionInstallation].Installation_INSPIRE_ID, [3_ProductionInstallation].mainActivityCode FROM "
                "3_ProductionInstallation;")
    with open(path_ied_match_facil_insp_inst_insp, 'w', newline='', encoding="utf-8") as f:
        writer = csv.writer(f, delimiter=";", quotechar='"')
        writer.writerow([column[0] for column in cur.description])
        for row in cur.fetchall():
            writer.writerow(row)
    cur.close()

    # Close MS ACCESS database connection
    conn.close()


if __name__ == "__main__":
    query_ied_data(
        path_ied_database=snakemake.input[0],
        path_ied_etsidentifiers=snakemake.output.ied_etsidentifiers,
        path_ied_general_facility_infos=snakemake.output.ied_general_facility_infos,
        path_ied_waste=snakemake.output.ied_waste,
        path_ied_pollutants=snakemake.output.ied_pollutants,
        path_ied_waste_water=snakemake.output.ied_waste_water,
        path_ied_match_facil_insp_inst_insp=snakemake.output.ied_match_facil_insp_inst_insp
    )
