localrules: eurostat_data_tsv, ch_data_xlsx

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
        src = script_dir + "eurostat/annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-nrg_bal_c.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = config["mapping-keys"]["eurostat"]["category-names"],
        carrier_names = config["mapping-keys"]["eurostat"]["carrier-names"]
    output: "build/data/annual-energy-balances.csv"
    params:
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"
