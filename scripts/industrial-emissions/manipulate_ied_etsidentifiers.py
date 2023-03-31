import pandas as pd
import re


def manipulate_ied_etsidentifiers(path_ied_etsidentifiers: str,
                                  path_ied_etsidentifiers_wconvetsid: str,
                                  reportingYear: int) -> None:

    """
    Background and Purpose of Script:
    For most installations the IED dataset reports ETSIdentifiers which are meant to establish a link between IED
    installations and ETS installations. Each ETS installation has a unique ID in the ISO format (Country code + 15
    digits, e.g. AT000000000000073) which is typically reported in the field ETSIdentifier. However, for some
    countries only the non-zero digits of the ISO code are reported in the field ETSIdentifiers (e.g., 73). This
    requires extending the code to the ISO format (AT000000000000073). And for other countries the ETSIdentifier does
    not report the installation's ISO ID, but its permitID, eprtrID or parts of either. Permit IDs and EPRTR IDs are
    national IDs that can have various formats. Moreover, over the years countries change which ETSIdentifier they
    provide (we observe that the ISO codes are increasingly used).
    This script converts the reported ETSidentifiers to a format that matches the IDs reported in the ETS dataset (ISO
    ID, permitID, eprtrID). The new, added column to the IED table "ETSIdentifiers" is called conv_ETSIdentifiers.
    Later, to establish a matching between IED and ETS installations, the IED installation's conv_ETSIdentifier will
    be compared to all entries in the ETS dataset columns ISO_id, permitID, eprtrID to identify potential matches.

    This script is based on the countries' reporting style of 2020 and needs to be redone, if IED data from a
    different year shall be used.

    A manual check of the datasets has identified the following reporting styles for the ETSIdentifier in 2020 (X
    represents a digit):

    AT provides ISO Codes
    BE provides ISO Codes for some installations, permitIDs (Format: BRXXXX) for other
    DE provides ISO Codes for some installtions (Brandenburg), permitIDs for other (Saarland, Format: XXXXX-XXX
        or XXXXX-XXXX)
    EE provides ISO Codes
    ES provides Codes (Format: 1-,2-,3-,4-,6-,8-digit numbers) that can be extended to ISO Codes for some
        installations, permitIDs (Format: ESXXXXXXXXXXXX) for others.
    FR provides ISO Codes for some installations, and for others Codes (Format: FR-new-XXXXXXXX) that can be reduced
        to permitIDs (Format: XXXXXXXX)
    GB provides permitIDs (Format: UK-E-XXXXX or UK-W-XXXXX or UK-N-XXXXX or UK-S-XXXXX)
    GR provides ISO Codes
    HR provides ISO Codes
    IE provides ISO Codes, however some codes have to be shortened (from Format: "IEXXXXXXXXXXXXXXX_PXXXX") and two
        installations contain two ISO codes (Format: "IEXXXXXXXXXXXXXXX , IEXXXXXXXXXXXXXXX"). Here, currently only
        the first is kept.
    LU provides ISO Codes
    PL provides Codes (Format: 1-,2-,3-,6-digit numbers) that can be extended to ISO Codes
    PT provides Codes (Format: 1-,2-,3-,6-digit numbers) that can be extended to ISO Codes for some installations,
        eprtrIDs (Format: XXXXXXXXX) for some installations, and permitIDs (Only 1 installation: 319.01III) for some
        installations
    SE provides ISO Codes

    Note: For the other countries no IED data exist or no ETSIdentifiers have been specified.
    Note: Some countries provide more installations with ETSIdentifiers, where the ETSIdentifier's format could not be
          resolved into either an ISO Code, permitID, or eprtrID.

    Parameters:
        path_ied_etsidentifiers (str):
            Input csv file with ETS Identifiers as reported in the IED dataset.
        path_ied_etsidentifiers_wconvetsid (str):
            Output csv file with ETS Identifiers transformed to a format that matches the ETS dataset's ISO ID, Permit
            ID, or E-PRTR ID.
        reportingYear (int):
            Choosing the reporting year for which the ETS Identifiers shall be transformed.

    """

    df = pd.read_csv(path_ied_etsidentifiers, sep=';', quotechar='"')
    df.insert(len(df.columns), "conv_ETSIdentifier", "")

    # Drop installation entries from all other years
    df.drop(df[df.reportingYear != reportingYear].index, inplace=True)

    for index, row in df.iterrows():

        if row["countryCode"] == "AT" or "EE" or "GB" or "GR" or "HR" or "LU" or "SE":
            df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]

        if row["countryCode"] == "BE":
            if bool(re.search("^BE\d{15,15}$", row["ETSIdentifier"])) or bool(re.search("^BR\d{6,6}$",
                                                                                        row["ETSIdentifier"])):
                df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]

        if row["countryCode"] == "DE":
            if bool(re.search("^DE\d{15,15}$", row["ETSIdentifier"])) or bool(re.search("^\d{5,5}-\d{3,4}$",
                                                                                        row["ETSIdentifier"])):
                df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]

        if row["countryCode"] == "ES":
            if row["ETSIdentifier"].isdigit():
                conv_to_ETSISO = "ES" + ("0" * (15 - len(row["ETSIdentifier"]))) + row["ETSIdentifier"]
                df.loc[index, "conv_ETSIdentifier"] = conv_to_ETSISO
            elif bool(re.search("^ES\d{12,12}$", row["ETSIdentifier"])):
                df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]

        if row["countryCode"] == "FR":
            if bool(re.search("FR-new-", row["ETSIdentifier"])):
                red_to_permitID = re.sub("\D", "", row["ETSIdentifier"])
                df.loc[index, "conv_ETSIdentifier"] = red_to_permitID
            elif bool(re.search("^FR\d{15,15}$", row["ETSIdentifier"])):
                df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]

        if row["countryCode"] == "IE":
            if bool(re.search("IE\d{15,15}", row["ETSIdentifier"])):
                list_of_ISOs = re.findall("IE\d{15,15}", row["ETSIdentifier"])
                df.loc[index, "conv_ETSIdentifier"] = list_of_ISOs[0]

        if row["countryCode"] == "PL":
            conv_to_ETSISO = "PL" + ("0" * (15 - len(row["ETSIdentifier"]))) + row["ETSIdentifier"]
            df.loc[index, "conv_ETSIdentifier"] = conv_to_ETSISO

        if row["countryCode"] == "PT":
            if bool(re.search("^\d{9,9}$", row["ETSIdentifier"])):
                df.loc[index, "conv_ETSIdentifier"] = row["ETSIdentifier"]
            if row["ETSIdentifier"].isdigit() and len(row["ETSIdentifier"]) <= 8:
                conv_to_ETSISO = "PT" + ("0" * (15 - len(row["ETSIdentifier"]))) + row["ETSIdentifier"]
                df.loc[index, "conv_ETSIdentifier"] = conv_to_ETSISO

    df.to_csv(path_ied_etsidentifiers_wconvetsid, index=False, sep=';', quotechar='"')


if __name__ == "__main__":
    manipulate_ied_etsidentifiers(
        path_ied_etsidentifiers=snakemake.input[0],
        path_ied_etsidentifiers_wconvetsid=snakemake.output[0],
        reportingYear=snakemake.params.reportingYear_IED
    )
