"""Rules to compute a unit's shares of the national electricity demand in each industrial sector."""


rule download_ets_data:
    message: "Download raw emissions as reported to the EU Emission Trading System."
    params: url = config["data-sources"]["emissions-ets"]
    output: protected("data/automatic/raw-emissions-ets.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_ets_data:
    message: "Unzip raw emissions as reported to the EU Emission Trading System."
    input: rules.download_ets_data.output
    output:
        ets_installation = "build/data/industrial-emissions/ets-data/installation.csv",
        ets_compliance = "build/data/industrial-emissions/ets-data/compliance.csv"
    shadow: "minimal"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} -d build/data/industrial-emissions/ets-data"


rule download_ied_data:
    message: "Download raw emissions as reported under the IED Directive."
    params: url = config["data-sources"]["emissions-ied"]
    output: protected("data/automatic/raw-emissions-ied.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_ied_data:
    message: "Unzip raw emissions as reported under the IED Directive."
    input: rules.download_ied_data.output
    output: "build/data/industrial-emissions/ied-data/1215_Public_Product_Full Access_draft_v19_April_2022_v6.accdb"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} -d build/data/industrial-emissions/ied-data"


rule query_ied_data:
    message: "Extract queries from the ied dataset (accdb format). This rule can only be executed on Windows machines."
    input: rules.unzip_ied_data.output
    output:
        ied_etsidentifiers = "build/data/industrial-emissions/ied-data/ied-etsidentifiers.csv",
        ied_general_facility_infos = "build/data/industrial-emissions/ied-data/ied-general-facility-infos.csv",
        ied_waste = "build/data/industrial-emissions/ied-data/ied-waste.csv",
        ied_pollutants = "build/data/industrial-emissions/ied-data/ied-pollutants.csv",
        ied_waste_water = "build/data/industrial-emissions/ied-data/ied-waste-water.csv",
        ied_match_facil_insp_inst_insp = "build/data/industrial-emissions/ied-data/ied-match-facil-insp-inst-insp.csv"
    conda: "../envs/query-msaccess.yaml"
    script: "../scripts/industrial-emissions/query_ied_data.py"


rule manipulate_ied_etsidentifiers:
    message:
        "Manipulate IED_ETSIdentifiers in preparation for matching of datasets. Script has to be adapted when using"
        "IED data from other years than 2020."
    input: rules.query_ied_data.output.ied_etsidentifiers # When importing these from a windows machine, make sure
        # that modification time stemps are more recent than those of unzip_ied_data's output
    params: reportingYear_IED = 2020 # ASSUME ratio emissions of sector j in unit i / national emissions of sector j
        # constant over years. Script has to be adapted when using IED data from other years than 2020.
    output: "build/data/industrial-emissions/ied-data/ied-etsidentifiers-converted-etsid.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/industrial-emissions/manipulate_ied_etsidentifiers.py"


rule preprocess_emission_data:
    message:
        "Preprocess ETS data. Preprocess IED data. Identify plants that occur in both datasets (matching). Join both"
        "datasets. Convert plants' activity codes to EUROSTAT format."
    input:
        ied_etsidentifiers = rules.manipulate_ied_etsidentifiers.output[0],
        ied_general_facility_infos = rules.query_ied_data.output.ied_general_facility_infos,
        ied_waste = rules.query_ied_data.output.ied_waste,
        ied_pollutants = rules.query_ied_data.output.ied_pollutants,
        ied_waste_water = rules.query_ied_data.output.ied_waste_water,
        ied_match_facil_insp_inst_insp = rules.query_ied_data.output.ied_match_facil_insp_inst_insp,
        ets_installation = rules.unzip_ets_data.output.ets_installation,
        ets_compliance = rules.unzip_ets_data.output.ets_compliance,
    params:
        reportingYear_ETS = 2020, # ASSUME ratio (emissions of sector j in unit i / national emissions of sector j)
            # constant over years. Reporting year for the IED dataset is set in rule manipulate_ied_etsidentifiers
        mapping_nace_ets = "config/mapping-nace-activityidets.csv",
        mapping_nace_eurostat = "config/mapping-nace-eurostat-categories.csv"
    output:
        industrial_emission_data_master = "data/industrial-emission-data-master.csv",
        industrial_emission_data_installation_details = "data/industrial-emission-data-installation-detail.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/industrial-emissions/preprocess_emission_data.py"
