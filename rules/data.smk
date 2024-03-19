"Rules downloading and preprocessing data that is used in different places of the workflow"
# TODO consider merging shape.smk into this

"Rules regarding CH Data:"

rule download_ch_energy_data:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"][f"swiss-{wildcards.dataset}"]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        dataset = "energy-balance|industry-energy-balance|end-use"
    localrule: True
    shell: "curl -sLo {output} {params.url}"


"Rules regarding Eurostat Data:"


rule download_eurostat_energy_data:
    message: "Download {wildcards.dataset} Eurostat data from euro-calliope datasets"
    params:
        url = lambda wildcards: config["data-sources"][f"eurostat-{wildcards.dataset}"]
    wildcard_constraints:
        dataset = "energy-balance|hh-end-use"
    conda: "../envs/shell.yaml"
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    localrule: True
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        energy_balance = "data/automatic/eurostat-energy-balance.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = "config/energy-balances/energy-balance-category-names.csv",
        carrier_names = "config/energy-balances/energy-balance-carrier-names.csv"
    output: "build/data/annual-energy-balances.csv"
    params:
        first_year = 2000
    conda: "../envs/default.yaml"
    script: "../scripts/data/annual_energy_balance.py"


"Rules regarding JRC-IDEES Data:"


rule download_jrc_idees_zipped:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: url = config["data-sources"]["jrc-idees"]
    output: protected("data/automatic/jrc-idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sLo {output} '{params.url}'"


rule jrc_idees_unzipped:
    message: "Unzip all JRC-IDEES {wildcards.sector} sector country data"
    input:
        "data/automatic/jrc-idees/{country_code}.zip"
    params: sector_title_case = lambda wildcards: wildcards.sector.title()
    wildcard_constraints:
        sector = "transport"
    output: temp("build/data/jrc-idees/{sector}/unprocessed/JRC-IDEES-2015_{sector}_{country_code}.xlsx")
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    localrule: True
    shell: "unzip -j {input} -d build/data/jrc-idees/{wildcards.sector}/unprocessed/"

rule jrc_idees_unzipped_heat:  # TODO merge this rule with the previous one. Careful to name confusion between 'heat' and 'Tertiary'
    message: "Unzip all JRC-IDEES heat sector country data"
    input:
        "data/automatic/jrc-idees/{country_code}.zip"
    output: temp("build/data/jrc-idees/heat/unprocessed/JRC-IDEES-2015_Tertiary_{country_code}.xlsx")
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    localrule: True
    shell: "unzip -j {input} -d build/data/jrc-idees/heat/unprocessed/"

"EU28 county codes used for downloading JRC-IDEES"
EU28 = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]

rule jrc_idees_transport_processed:
    message: "Process {wildcards.dataset} transport data from JRC-IDEES to be used in understanding current and future transport demand"
    input:
        data = expand(
            "build/data/jrc-idees/transport/unprocessed/JRC-IDEES-2015_Transport_{country_code}.xlsx",
            country_code=EU28
        )
    output: "build/data/jrc-idees/transport/processed-{dataset}.csv"
    params:
        vehicle_type_names = config["parameters"]["transport"]["vehicle-type-names"],
    wildcard_constraints:
        dataset = "road-energy|road-distance|road-vehicles"
    conda: "../envs/default.yaml"
    script: "../scripts/transport/jrc_idees.py"

rule jrc_idees_heat_processed: # TODO merge this rule with the previous one. Careful to name confusion between 'heat' and 'Tertiary'
    message: "Process tertiary heat data from JRC-IDEES"
    input:
        data = expand(
            "build/data/jrc-idees/heat/unprocessed/JRC-IDEES-2015_Tertiary_{country_code}.xlsx",
            country_code=EU28
        )
    output: "build/data/jrc-idees/heat/commercial/jrc_idees_processed_energy.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/heat/jrc_idees.py"