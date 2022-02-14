"""Rules to generate electricity demand time series."""

localrules: download_raw_load


rule download_raw_load:
    message: "Download raw load."
    params: url = config["data-sources"]["load"]
    output: protected("data/automatic/raw-load-data.csv")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        script = script_dir + "demand/national_load.py",
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = config["quality-control"]["load"],
        countries = config["scope"]["spatial"]["countries"]
    output: "build/data/electricity-demand-national.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/demand/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        script = script_dir + "demand/load.py",
        units = rules.units.output[0],
        demand_per_unit = rules.potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/timeseries/demand/electricity.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/load.py"


rule electricity_demand_techs_and_locations_template:
    message: "Create {wildcards.resolution} tech definition file from template."
    input:
        script = script_dir + "template_techs.py",
        template = techs_template_dir + "demand/electricity.yaml",
        locations = rules.locations_template.output.csv,
        timeseries_data = rules.electricity_load.output
    params:
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/demand/electricity.yaml"
    script: "../scripts/template_techs.py"
