localrules: ch_data_xlsx


rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][wildcards.dataset]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


rule ch_annual_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        script = script_dir + "ch-stats/annual_energy_balance.py",
        ch_energy_balance_excel = "data/automatic/ch-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/annual-energy-balances.nc"
    script: "../scripts/ch-stats/annual_energy_balance.py"


rule ch_annual_industry_subsector_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        script = script_dir + "ch-stats/industry_subsectors.py",
        ch_industry_excel = "data/automatic/ch-industry-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/industry-energy-balance.nc"
    script: "../scripts/ch-stats/industry_subsectors.py"


rule ch_annual_building_end_use_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        script = script_dir + "ch-stats/building_heat.py",
        ch_end_use_excel = "data/automatic/ch-end-use.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/building-heat-energy-balance.nc"
    script: "../scripts/ch-stats/building_heat.py"


rule ch_gross_added_value:
    message: "Process Swiss gross added value spreadsheets"
    input:
        script = script_dir + "ch-stats/gross_added_value.py",
        ch_gva_excel = "data/automatic/ch-gva.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/gross_added_value.csv"
    script: "../scripts/ch-stats/gross_added_value.py"


rule ch_annual_transport_energy_balance:
    message: "Process Swiss energy balance spreadsheets for transport data"
    input:
        script = script_dir + "ch-stats/transport.py",
        ch_energy_balance_excel = "data/automatic/ch-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/transport-energy-balance.nc"
    script: "../scripts/ch-stats/transport.py"
