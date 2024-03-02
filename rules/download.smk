localrules: ch_data_xlsx, eurostat_annual_energy_balances, jrc_idees_zipped


"Rules regarding CH Data:"
rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][wildcards.dataset]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


"Rules regarding Eurostat Data:"

rule eurostat_annual_energy_balances:
    message: "Get Annual Energy Balances from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-energy-balance"]
    output: protected("data/automatic/eurostat-energy-balance.tsv.gz")
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
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"


"Rules regarding JRC-IDEES Data:"

"EU28 county codes used for downloading JRC-IDEES"
EU28 = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]

rule jrc_idees_zipped:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: url = config["data-sources"]["jrc-idees"]
    output: protected("data/automatic/jrc-idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule jrc_idees_unzipped:
    message: "Unzip all JRC-IDEES {wildcards.sector} sector country data"
    input:
        countries = [
            f"data/automatic/jrc-idees/{country_code}.zip"
            for country_code in [
                pycountry.countries.lookup(country).alpha_2 if pycountry.countries.lookup(country).name != "Greece" else "EL"
                for country in config["scope"]["spatial"]["countries"]
            ]
            if country_code in EU28
        ]
    params: sector_title_case = lambda wildcards: wildcards.sector.title()
    wildcard_constraints:
        sector = "((industry)|(transport)|(tertiary))"
    output: temp(directory("build/data/jrc-idees/{sector}/unprocessed"))
    conda: "../envs/shell.yaml"
    shell: "unzip 'data/automatic/jrc-idees/*.zip' '*{params.sector_title_case}*' -d {output}"


rule jrc_idees_transport_processed:
    message: "Process {wildcards.dataset} transport data from JRC-IDEES to be used in understanding current and future transport demand"
    input:
        unprocessed_data = "build/data/jrc-idees/transport/unprocessed"
    output: "build/data/jrc-idees/transport/processed-{dataset}.csv"
    wildcard_constraints:
        dataset = "((road-energy)|(road-distance)|(road-vehicles))"
    conda: "../envs/default.yaml"
    script: "../scripts/transport/jrc_idees.py"
