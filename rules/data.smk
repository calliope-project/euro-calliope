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
    shell: "curl -sSLo {output} {params.url}"


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
    shell: "curl -sSLo {output} {params.url}"


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
    shell: "curl -sSLo {output} '{params.url}'"


def jrc_to_euro_calliope_sector(sector: str):
    if sector == "transport":
        return "Transport"
    elif sector == "heat":
        return "Tertiary"
    else:
        raise ValueError(f"Unknown sector {sector}.")


rule jrc_idees_unzipped:
    message: "Unzip all JRC-IDEES {wildcards.sector} sector country data"
    input:
        "data/automatic/jrc-idees/{country_code}.zip"
    params:
        file_name = lambda wildcards: f"JRC-IDEES-2015_{jrc_to_euro_calliope_sector(wildcards.sector)}_{wildcards.country_code}.xlsx"
    wildcard_constraints:
        sector = "transport|heat"
    output: temp("build/data/jrc-idees/{sector}/unprocessed/{country_code}.xlsx")
    conda: "../envs/shell.yaml"
    shadow: "minimal"
    localrule: True
    shell: """
    unzip -j {input} -d build/data/jrc-idees/{wildcards.sector}/unprocessed/
    mv build/data/jrc-idees/{wildcards.sector}/unprocessed/{params.file_name} {output}
    """

"EU28 county codes used for downloading JRC-IDEES"
JRC_IDEES_SCOPE = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]

rule download_gridded_temperature_data:
    message: "Download gridded temperature data"
    params: url = config["data-sources"]["gridded-temperature-data"]
    output: protected("data/automatic/gridded-weather/temperature.nc")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"
