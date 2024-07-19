"""Rules to generate electricity demand time series."""


rule download_raw_load:
    message: "Download raw load."
    params: url = config["data-sources"]["load"]
    output: protected("data/automatic/raw-load-data.csv")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sLo {output} '{params.url}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = config["quality-control"]["load"],
        countries = config["scope"]["spatial"]["countries"]
    output:
        csv = "build/data/electricity-demand-national.csv"
    log: "logs/electricity_load_national.log"
    conda: "../envs/default.yaml"
    script: "../scripts/demand/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        units = rules.units.output[0],
        demand_per_unit = rules.potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/timeseries/demand/electricity.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/load.py"
