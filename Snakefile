PYTHON = "PYTHONPATH=./ python"

URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"

TECHNICAL_POTENTIAL_COUNTRIES = "src/data/national-technical-potential.geojson"

rule all:
    message: "Generate Euro Calliope."
    input:
        "model/model.done"


rule countries:
    message: "Generate locations for all countries."
    input: TECHNICAL_POTENTIAL_COUNTRIES
    output: "model/countries.yaml"
    conda: "src/envs/geo.yaml"
    script: "src/locations.py"


rule raw_load:
    message: "Download raw load."
    output: protected("src/data/automatic/raw-load-data.csv")
    shell: "curl -sLo {output} '{URL_LOAD}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input: rules.raw_load.output
    output: "src/data/generated/electricity-demand-national.csv"
    params:
        number_rows_valid = 10654293, # see https://github.com/Open-Power-System-Data/time_series/issues/22
        year = 2017
    script: "src/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location."
    input:
        units = TECHNICAL_POTENTIAL_COUNTRIES,
        national_load = rules.electricity_load_national.output[0]
    output: "model/electricity-demand.csv"
    conda: "src/envs/geo.yaml"
    script: "src/load.py"


rule model:
    message: "Generate Euro Calliope"
    input:
        rules.countries.output,
        rules.electricity_load.output
    output: touch("model/model.done")


rule clean: # removes all generated results
    shell:
        """
        rm -r model/
        rm -r src/data/generated/
        """


rule test:
    message: "Run tests"
    input:
        rules.model.output
    shell:
        "py.test"
