URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"

LOCATIONS = "src/data/national-technical-potential.geojson"
EEZ = "src/data/eez-in-europe.geojson"
SHARED_COAST = "src/data/national-shared-coast.csv"
LAND_ELIGIBILITY = "src/data/national-eligibility.csv"


rule all:
    message: "Generate Euro Calliope and run tests."
    input:
        "model/model.done",
        "model/test-report.html"


rule locations:
    message: "Generate locations."
    input:
        land_eligibility_km2 = LAND_ELIGIBILITY
    output: "model/locations.yaml"
    conda: "src/envs/default.yaml"
    script: "src/locations.py"


rule capacity_factors:
    message: "Generate capacityfactor time series disaggregated by location for {wildcards.technology}."
    input:
        locations = LOCATIONS,
        ids = "src/data/capacityfactors/{technology}-ids.tif",
        timeseries = "src/data/capacityfactors/{technology}-timeseries.nc"
    output: "model/capacityfactors-{technology,^((?!offshore).)*}.csv"
    conda: "src/envs/geo.yaml"
    script: "src/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location for wind-offshore."
    input:
        eez = EEZ,
        shared_coast = SHARED_COAST,
        ids = "src/data/capacityfactors/wind-offshore-ids.tif",
        timeseries = "src/data/capacityfactors/wind-offshore-timeseries.nc"
    output: "model/capacityfactors-wind-offshore.csv"
    conda: "src/envs/geo.yaml"
    script: "src/capacityfactors_offshore.py"


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
        year = 2016
    conda: "src/envs/default.yaml"
    script: "src/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location."
    input:
        units = LOCATIONS,
        national_load = rules.electricity_load_national.output[0]
    output: "model/electricity-demand.csv"
    conda: "src/envs/geo.yaml"
    script: "src/load.py"


rule model:
    message: "Generate Euro Calliope"
    input:
        rules.locations.output,
        rules.electricity_load.output,
        expand(
            "model/capacityfactors-{technology}.csv",
            technology=["rooftop-pv", "open-field-pv", "wind-onshore", "wind-offshore"]
        )
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
    output: "model/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
