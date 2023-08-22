localrules: ch_data_xlsx

rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][wildcards.dataset]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"

rule ch_annual_transport_energy_balance:
    message: "Process Swiss energy balance spreadsheets for transport data"
    input:
        script = script_dir + "ch-stats/transport.py",
        ch_energy_balance_excel = "data/automatic/ch-energy-balance.xlsx"
    conda: "../envs/default.yaml"
    output: "build/data/ch-stats/transport-energy-balance.nc"
    script: "../scripts/ch-stats/transport.py"
