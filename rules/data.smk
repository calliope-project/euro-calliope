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


rule download_raw_population_zipped:
    message: "Download population data."
    output:
        protected("data/automatic/raw-population-data.zip")
    params: url = config["data-sources"]["population"]
    conda: "../envs/shell.yaml"
    shell: "curl -sSLo {output} '{params.url}'"


rule raw_population_unzipped:
    message: "Extract population data TIF."
    input: rules.download_raw_population_zipped.output
    output: temp("build/JRC_1K_POP_2018.tif")
    conda: "../envs/shell.yaml"
    shell: "unzip {input} '*.tif' -d ./build/"
