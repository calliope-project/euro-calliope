"""Rules to process input data."""

configfile: "./config/default.yaml"
localrules: eurostat_data_tsv, ch_data_xlsx
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule eurostat_data_tsv:
    message: "Get various datasets from Eurostat"
    output:
        energy_balance = protected("data/automatic/annual_energy_balances.tsv.gz"),
        hh_end_use = protected("data/automatic/hh_end_use.tsv.gz"),
        freight = protected("data/automatic/freight.tsv.gz"),
        employees = protected("data/automatic/employees.tsv.gz"),
        gva = protected("data/automatic/gva.tsv.gz"),
        dwellings = protected("data/automatic/dwellings.tsv.gz"),
    shell:
        """
        curl -sLo {output.energy_balance} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][energy-balance]}.tsv.gz'
        curl -sLo {output.hh_end_use} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][hh-end-use]}.tsv.gz'
        curl -sLo {output.freight} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][freight]}.tsv.gz'
        curl -sLo {output.employees} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][employees]}.tsv.gz'
        curl -sLo {output.gva} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][gva]}.tsv.gz'
        curl -sLo {output.dwellings} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][dwellings]}.tsv.gz'
        """


rule ch_data_xlsx:
    message: "Get Swiss annual energy balances and household end uses"
    output:
        energy_balance = protected("data/automatic/ch_annual_energy_balances.xlsx"),
        industry_energy_balance = protected("data/automatic/ch_annual_industry_energy_balances.xlsx"),
        end_use = protected("data/automatic/ch_hh_end_use.xlsx"),
        gva = protected("data/automatic/ch_gva.xlsx")
    shell:
        """
        curl -sLo {output.energy_balance} '{config[data-sources][swiss-stat][energy-balance]}'
        curl -sLo {output.industry_energy_balance} '{config[data-sources][swiss-stat][industry-energy-balance]}'
        curl -sLo {output.end_use} '{config[data-sources][swiss-stat][end-use]}'
        curl -sLo {output.gva} '{config[data-sources][swiss-stat][gva]}'
        """


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat and Switzerland-specific data"
    input:
        src = script_dir + "annual_energy_balance.py",
        energy_balance = rules.eurostat_data_tsv.output.energy_balance,
        ch_energy_balance = rules.ch_data_xlsx.output.energy_balance,
        ch_industry_energy_balance = rules.ch_data_xlsx.output.industry_energy_balance,
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    shadow: "minimal"
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
