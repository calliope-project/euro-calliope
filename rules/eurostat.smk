localrules: eurostat_data_tsv

rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-base-url"].format(dataset=wildcards.dataset)
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat data"
    input:
        src = script_dir + "eurostat/annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-nrg_bal_c.tsv.gz",
    output: "build/data/eurostat/annual-energy-balances.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/annual_energy_balance.py"


rule household_building_heat:
    message: "Process household building heat end use annual energy balances from Eurostat"
    input:
        src = script_dir + "eurostat/household_building_heat.py",
        household_end_use_energy_balance = "data/automatic/eurostat-nrg_d_hhq.tsv.gz",
    output: "build/data/eurostat/household-building-heat-end-use-energy-balances.nc"
    conda: "../envs/default.yaml"
    script: "../scripts/eurostat/household_building_heat.py"


rule sectoral_annual_energy_balances:
    message: "Blend swiss and Eurostat energy balances for the {wildcards.sector} sector"
    input:
        script = script_dir + "eurostat/blend_and_rename_per_sector.py",
        eurostat_energy_balances="build/data/eurostat/annual-energy-balances.nc",
        ch_energy_balances="build/data/ch-stats/annual-energy-balances.nc"
    params:
        carrier_name_mapping = config["statistical-code-mapping"]["eurostat-to-sectoral-carrier-names"],
        category_name_mapping = config["statistical-code-mapping"]["eurostat-to-sectors"]
    wildcard_constraints:
        sector = "commercial|household"
    conda: "../envs/default.yaml"
    output: "build/data/eurostat/annual-{sector}-energy-balances.nc"
    script: "../scripts/eurostat/blend_and_rename_per_sector.py"


use rule sectoral_annual_energy_balances as transport_subsector_annual_energy_balances with:
    input:
        script = script_dir + "eurostat/blend_and_rename_per_sector.py",
        eurostat_energy_balances="build/data/eurostat/annual-energy-balances.nc",
        ch_energy_balances="build/data/ch-stats/transport-energy-balance.nc"
    wildcard_constraints:
        sector = "road-transport|aviation|marine|rail"


use rule sectoral_annual_energy_balances as industry_subsector_annual_energy_balances with:
    input:
        script = script_dir + "eurostat/blend_and_rename_per_sector.py",
        eurostat_energy_balances="build/data/eurostat/annual-energy-balances.nc",
        ch_energy_balances="build/data/ch-stats/industry-energy-balance.nc"
    params:
        carrier_name_mapping = config["statistical-code-mapping"]["eurostat-to-sectoral-carrier-names"],
        category_name_mapping = config["statistical-code-mapping"]["industry-eurostat-to-jrc-idees"]
    wildcard_constraints:
        sector = "industry"


