URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"

LOCATIONS = "data/national-technical-potential.geojson"
EEZ = "data/eez-in-europe.geojson"
SHARED_COAST = "data/national-shared-coast.csv"
LAND_ELIGIBILITY = "data/national-eligibility.csv"

localrules: all, raw_load, model, clean, copy_template


onstart:
    shell("mkdir -p build/logs")


rule all:
    message: "Generate Euro Calliope and run tests."
    input:
        "build/logs/model.done",
        "build/logs/test-report.html"


rule copy_template:
    message: "Copy file {wildcards.definition_file}.yaml from templates."
    input: "src/template/{definition_file}.yaml",
    output: "build/model/{definition_file}.yaml"
    shell: "cp {input} {output}"


rule locations:
    message: "Generate locations."
    input:
        src = "src/locations.py",
        land_eligibility_km2 = LAND_ELIGIBILITY
    output: "build/model/locations.yaml"
    conda: "envs/default.yaml"
    script: "src/locations.py"


rule capacity_factors:
    message: "Generate capacityfactor time series disaggregated by location for {wildcards.technology}."
    input:
        src = "src/capacityfactors.py",
        locations = LOCATIONS,
        ids = "data/capacityfactors/{technology}-ids.tif",
        timeseries = "data/capacityfactors/{technology}-timeseries.nc"
    wildcard_constraints:
        technology = "((wind-onshore)|(rooftop-pv)|(open-field-pv))"
    output: "build/model/capacityfactors-{technology}.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location for wind-offshore."
    input:
        src = "src/capacityfactors_offshore.py",
        eez = EEZ,
        shared_coast = SHARED_COAST,
        ids = "data/capacityfactors/wind-offshore-ids.tif",
        timeseries = "data/capacityfactors/wind-offshore-timeseries.nc"
    output: "build/model/capacityfactors-wind-offshore.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors_offshore.py"


rule raw_load:
    message: "Download raw load."
    output: protected("data/automatic/raw-load-data.csv")
    shell: "curl -sLo {output} '{URL_LOAD}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        src = "src/national_load.py",
        load = rules.raw_load.output
    output: "build/data/electricity-demand-national.csv"
    params:
        number_rows_valid = 10654293, # see https://github.com/Open-Power-System-Data/time_series/issues/22
        year = 2016
    conda: "envs/default.yaml"
    script: "src/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location."
    input:
        src = "src/load.py",
        units = LOCATIONS,
        national_load = rules.electricity_load_national.output[0]
    output: "build/model/electricity-demand.csv"
    conda: "envs/geo.yaml"
    script: "src/load.py"


rule model:
    message: "Generate Euro Calliope"
    input:
        "build/model/interest-rate.yaml",
        "build/model/link-techs.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        rules.locations.output,
        rules.electricity_load.output,
        expand(
            "build/model/capacityfactors-{technology}.csv",
            technology=["rooftop-pv", "open-field-pv", "wind-onshore", "wind-offshore"]
        )
    output: touch("build/logs/model.done")


rule clean: # removes all generated results
    shell:
        """
        rm -r build/
        """


rule test:
    message: "Run tests"
    input:
        rules.model.output
    conda: "envs/test.yaml"
    output: "build/logs/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
