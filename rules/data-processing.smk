"""Rules to process input data."""

configfile: "./config/default.yaml"
localrules: eurostat_data_tsv, ch_data_xlsx
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-base-url"].format(dataset=wildcards.dataset)
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.url}"


rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][wildcards.dataset]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat and Switzerland-specific data"
    input:
        src = script_dir + "annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-nrg_bal_c.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = config["data-sources"]["energy-balance-category-names"],
        carrier_names = config["data-sources"]["energy-balance-carrier-names"]
    output: "build/data/annual-energy-balances.csv"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/annual_energy_balance.py"

rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        src = "scripts/annual_heat_demand.py",
        hh_end_use = rules.eurostat_data_tsv.output.hh_end_use,
        ch_end_use = rules.ch_data_xlsx.output.end_use,
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "data/commercial/jrc_idees_processed_energy.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    params:
        countries = config["scope"]["countries"],
        heat_tech_params = config["parameters"]["heat-end-use"]
    conda: "../envs/default.yaml"
    output:
        demand=temp("build/annual_heat_demand.csv"),
        electricity=temp("build/annual_heat_electricity_consumption.csv"),
    script: "../scripts/annual_heat_demand.py"
