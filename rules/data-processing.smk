"""Rules to process input data."""

configfile: "./config/default.yaml"
localrules: eurostat_data_tsv, ch_data_xlsx
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        base_url = config["data-sources"]["eurostat-base-url"],
    output: protected("data/automatic/eurostat_{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.base_url}{wildcards.dataset}.tsv.gz"


rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][f"{wildcards.dataset}"],
    output: protected("data/automatic/ch_{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat and Switzerland-specific data"
    input:
        src = script_dir + "annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat_nrg_bal_c.tsv.gz",
        ch_energy_balance = "data/automatic/ch_energy_balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch_industry_energy_balance.xlsx",
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output: "`build/annual_energy_balances.csv`"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    shadow: "minimal"
    script: "../scripts/annual_energy_balance.py"
