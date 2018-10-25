PYTHON = "PYTHONPATH=./ python"

URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"

COUNTRIES = "src/data/national-technical-potential.geojson"
LAND_ELIGIBILITY = "src/data/national-eligibility.csv"


rule all:
    message: "Generate Euro Calliope and run tests."
    input:
        "model/model.done",
        "model/tests.done"


rule countries:
    message: "Generate locations for all countries."
    input:
        ids = COUNTRIES,
        land_eligibility = LAND_ELIGIBILITY
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
    conda: "src/envs/default.yaml"
    script: "src/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location."
    input:
        units = COUNTRIES,
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
    conda: "src/envs/test.yaml"
    output: touch("model/tests.done")
    shell:
        "py.test"
