URL_LOAD = "https://data.open-power-system-data.org/time_series/2018-06-30/time_series_60min_stacked.csv"

LOCATIONS = "data/{resolution}/units.geojson"
EEZ = "data/eez-in-europe.geojson"
SHARED_COAST = "data/{resolution}/shared-coast.csv"
LAND_ELIGIBILITY = "data/{resolution}/eligibility.csv"

localrules: all, raw_load, model, clean, scale_template
configfile: "config/default.yaml"


onstart:
    shell("mkdir -p build/logs")


rule all:
    message: "Generate Euro Calliope and run tests."
    input:
        "build/logs/national-model.done",
        "build/logs/regional-model.done",
        "build/logs/test-report.html"


rule scale_template:
    message: "Apply scaling factors to file {wildcards.definition_file}.yaml from templates."
    input:
        src = "src/scale_templates.py",
        template = "src/template/{definition_file}.yaml"
    params: scaling_factors = config["scaling-factors"]
    output: "build/model/{definition_file}.yaml"
    conda: "envs/default.yaml"
    script: "src/scale_templates.py"


rule locations:
    message: "Generate locations for {wildcards.resolution} resolution."
    input:
        src = "src/locations.py",
        shapes = LOCATIONS,
        land_eligibility_km2 = LAND_ELIGIBILITY
    params: scaling_factors = config["scaling-factors"]
    output: "build/model/{resolution}/locations.yaml"
    conda: "envs/geo.yaml"
    script: "src/locations.py"


rule capacity_factors:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        src = "src/capacityfactors.py",
        locations = LOCATIONS,
        ids = "data/capacityfactors/{technology}-ids.tif",
        timeseries = "data/capacityfactors/{technology}-timeseries.nc"
    wildcard_constraints:
        technology = "((wind-onshore)|(rooftop-pv)|(open-field-pv))"
    output: "build/model/{resolution}/capacityfactors-{technology}.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for wind-offshore."
    input:
        src = "src/capacityfactors_offshore.py",
        eez = EEZ,
        shared_coast = SHARED_COAST,
        ids = "data/capacityfactors/wind-offshore-ids.tif",
        timeseries = "data/capacityfactors/wind-offshore-timeseries.nc"
    output: "build/model/{resolution}/capacityfactors-wind-offshore.csv"
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
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        src = "src/load.py",
        units = LOCATIONS,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/model/{resolution}/electricity-demand.csv"
    conda: "envs/geo.yaml"
    script: "src/load.py"


rule link_neighbours:
    message: "Create links between all direct neighbours on {wildcards.resolution} resolution."
    input:
        src = "src/link_neighbours.py",
        units = LOCATIONS
    output: "build/model/{resolution}/link-all-neighbours.yaml"
    conda: "envs/geo.yaml"
    script: "src/link_neighbours.py"


rule model:
    message: "Generate Euro Calliope with {wildcards.resolution} resolution."
    input:
        "build/model/interest-rate.yaml",
        "build/model/link-techs.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        rules.locations.output,
        rules.electricity_load.output,
        rules.link_neighbours.output,
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=["rooftop-pv", "open-field-pv", "wind-onshore", "wind-offshore"]
        )
    output: touch("build/logs/{resolution}-model.done")


rule clean: # removes all generated results
    shell:
        """
        rm -r build/
        """


rule test:
    message: "Run tests for national and regional models."
    input:
        "build/logs/national-model.done",
        "build/logs/regional-model.done",
    conda: "envs/test.yaml"
    output: "build/logs/test-report.html"
    shell:
        "py.test --html={output} --self-contained-html"
