localrules: eurostat_annual_energy_balances

rule eurostat_annual_energy_balances:
    message: "Get Annual Energy Balances from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-energy-balance"]
    output: protected("data/automatic/eurostat-energy-balance.tsv.gz")
    shell: "curl -sLo {output} {params.url}"

rule annual_energy_balances:
    message: "Get annual energy balances from Eurostat"
    input:
        src = script_dir + "eurostat/annual_energy_balance.py",
        energy_balance = "data/automatic/eurostat-energy-balance.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = "config/energy-balances/energy-balance-category-names.csv",
        carrier_names = "config/energy-balances/energy-balance-carrier-names.csv"
    output: "build/annual-energy-balances.csv"
    params:
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"
