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
        cat_names = "config/energy-balances/energy_balance_category_names.csv",
        carrier_names = "config/energy-balances/energy_balance_carrier_names.csv"
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"


#this was the previous version
rule annual_energy_balances_old:
    message: "Process annual energy balances from Eurostat data"
    input:
        src = script_dir + "eurostat/annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-energy-balance.tsv.gz",
    output: "build/data/eurostat/annual-energy-balances.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"

rule transport_subsector_annual_energy_balances:
    message: "Blend swiss and Eurostat energy balances for the {wildcards.sector} sector"
    input:
        script = script_dir + "eurostat/blend_and_rename_per_sector.py",
        eurostat_energy_balances="build/data/eurostat/annual-energy-balances.nc",
        ch_energy_balances="build/data/ch-stats/transport-energy-balance.nc"
    params:
        carrier_name_mapping = config["statistical-code-mapping"]["eurostat-to-sectoral-carrier-names"],
        category_name_mapping = config["statistical-code-mapping"]["eurostat-to-sectors"]
    wildcard_constraints:
        sector = "road-transport"
    conda: "../envs/default.yaml"
    output: "build/data/eurostat/annual-{sector}-energy-balances.nc"
    script: "../scripts/eurostat/blend_and_rename_per_sector.py"



