"Rules downloading and processing Swiss statistics data"

rule download_ch_energy_data:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"][f"swiss-{wildcards.dataset}"]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    conda: "../envs/shell.yaml"
    wildcard_constraints:
        dataset = "energy-balance|industry-energy-balance|end-use"
    localrule: True
    shell: "curl -sLo {output} {params.url}"


rule ch_annual_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        ch_energy_balance_excel = "data/automatic/ch-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/annual-energy-balances.nc"
    script: "../scripts/ch-stats/annual_energy_balance.py"


rule ch_annual_industry_subsector_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        ch_industry_excel = "data/automatic/ch-industry-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/industry-energy-balance.nc"
    script: "../scripts/ch-stats/industry_subsectors.py"


rule ch_annual_building_end_use_energy_balance:
    message: "Process Swiss energy balance spreadsheets"
    input:
        ch_end_use_excel = "data/automatic/ch-end-use.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/building-heat-energy-balance.nc"
    script: "../scripts/ch-stats/building_heat.py"


rule ch_annual_transport_energy_balance:
    message: "Process Swiss energy balance spreadsheets for transport data"
    input:
        ch_energy_balance_excel = "data/automatic/ch-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/transport-energy-balance.nc"
    script: "../scripts/ch-stats/transport.py"
