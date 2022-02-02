import glob
from pathlib import Path

from snakemake.utils import validate

ALL_WIND_AND_SOLAR_TECHNOLOGIES = [
    "wind-onshore", "wind-offshore", "open-field-pv",
    "rooftop-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat"
]

configfile: "config/default.yaml"
validate(config, "config/schema.yaml")

include: "./rules/shapes.smk"
include: "./rules/wind-and-solar.smk"
include: "./rules/biofuels.smk"
include: "./rules/hydro.smk"
include: "./rules/sync.smk"
include: "./rules/model_template.smk"
localrules: all, download_raw_load, model, clean
localrules: download_entsoe_tyndp_zip
wildcard_constraints:
        resolution = "continental|national|regional"

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
__version__ = open(f"{root_dir}VERSION").readlines()[0].strip()
script_dir = f"{root_dir}scripts/"
test_dir = f"{root_dir}tests/"
model_test_dir = f"{test_dir}model"

onstart:
    shell("mkdir -p build/logs")
onsuccess:
     if "email" in config.keys():
         shell("echo "" | mail -s 'euro-calliope succeeded' {config[email]}")
onerror:
     if "email" in config.keys():
         shell("echo "" | mail -s 'euro-calliope failed' {config[email]}")


rule all:
    message: "Generate euro-calliope pre-built models and run tests."
    input:
        "build/logs/continental/test-report.html",
        "build/logs/national/test-report.html",
        "build/logs/continental/model.done",
        "build/logs/national/model.done",
        "build/logs/regional/model.done",
        "build/models/build-metadata.yaml"


rule all_tests:
    message: "Generate euro-calliope pre-built models and run all tests."
    input:
        "build/logs/continental/model.done",
        "build/logs/national/model.done",
        "build/logs/regional/model.done",
        "build/logs/continental/test-report.html",
        "build/logs/national/test-report.html",
        "build/logs/regional/test-report.html",
        "build/models/build-metadata.yaml"



rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        script = script_dir + "hydro_capacities.py",
        locations = rules.units.output[0],
        plants = rules.preprocess_hydro_stations.output[0]
    output: "build/data/{resolution}/hydro-capacities-mw.csv"
    conda: "envs/geo.yaml"
    script: "scripts/hydro_capacities.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        script = script_dir + "capacityfactors_hydro.py",
        capacities = rules.hydro_capacities.output[0],
        stations = "build/data/hydro-electricity-with-energy-inflow-{first_year}-{final_year}.nc".format(
            first_year = config["scope"]["temporal"]["first-year"],
            final_year = config["scope"]["temporal"]["final-year"]
        ),
        locations = rules.units.output[0]
    params:
        threshold = config["capacity-factors"]["min"]
    output:
        ror = "build/models/{resolution}/timeseries/supply/capacityfactors-hydro-ror.csv",
        reservoir = "build/models/{resolution}/timeseries/supply/capacityfactors-hydro-reservoir-inflow.csv"
    conda: "envs/geo.yaml"
    script: "scripts/capacityfactors_hydro.py"


rule download_raw_load:
    message: "Download raw load."
    params: url = config["data-sources"]["load"]
    output: protected("data/automatic/raw-load-data.csv")
    conda: "envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        script = script_dir + "national_load.py",
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = config["quality-control"]["load"],
        countries = config["scope"]["spatial"]["countries"]
    output: "build/data/electricity-demand-national.csv"
    conda: "envs/default.yaml"
    script: "scripts/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        script = script_dir + "load.py",
        units = rules.units.output[0],
        demand_per_unit = rules.potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/timeseries/demand/electricity-demand.csv"
    conda: "envs/geo.yaml"
    script: "scripts/load.py"


rule download_entsoe_tyndp_zip:
    message: "Download ENTSO-E ten-year network development plan (TYNDP) 2020 scenario dataset"
    params: url = config["data-sources"]["entsoe-tyndp"]
    output: protected("data/automatic/raw-entsoe-tyndp.xlsx.zip")
    conda: "envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule entsoe_tyndp_xlsx:
    message: "Unzip ENTSO-E TYNDP 2020 scenario dataset."
    input: rules.download_entsoe_tyndp_zip.output[0]
    shadow: "minimal"
    output: "build/data/national/TYNDP-2020-Scenario-Datafile.xlsx",
    conda: "envs/shell.yaml"
    shell: "unzip -o {input} 'TYNDP-2020-Scenario-Datafile.xlsx' -d build/data/national"


rule build_metadata:
    message: "Generate build metadata."
    input:
        script_dir + "metadata.py",
        "build/logs/continental/model.done",
        "build/logs/national/model.done",
        "build/logs/regional/model.done",
    params:
        config = config,
        version = __version__
    output: "build/models/build-metadata.yaml"
    conda: "envs/default.yaml"
    script: "scripts/metadata.py"


rule model:
    message: "Generate euro-calliope with {wildcards.resolution} resolution."
    input:
        "build/models/{resolution}/example-model.yaml",
        rules.electricity_load.output,
        rules.capacity_factors_hydro.output,
        rules.hydro_capacities.output,
        expand(
            "build/models/{{resolution}}/timeseries/supply/capacityfactors-{technology}.csv",
            technology=ALL_WIND_AND_SOLAR_TECHNOLOGIES
        ),
    output:
        log = "build/logs/{resolution}/model.done",
    shell:
        """
        touch {output.log}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r build/
        echo "Data downloaded to data/automatic/ has not been cleaned."
        """


rule test:
    message: "Run tests"
    input:
        "build/logs/{resolution}/model.done",
        test_dir = model_test_dir,
        tests = map(str, Path(model_test_dir).glob("**/test_*.py")),
        example_model = "build/models/{resolution}/example-model.yaml",
        capacity_factor_timeseries = expand(
            "build/models/{{resolution}}/timeseries/supply/capacityfactors-{technology}.csv",
            technology=ALL_WIND_AND_SOLAR_TECHNOLOGIES + ["hydro-ror", "hydro-reservoir-inflow"]
        )
    params:
        config = config,
        scenarios = lambda wildcards: config["parameters"]["test"]["scenarios"][wildcards.resolution],
        subset_time = lambda wildcards: config["parameters"]["test"]["subset_time"][wildcards.resolution],
    output: "build/logs/{resolution}/test-report.html"
    conda: "./envs/test.yaml"
    script: "./tests/model/test_runner.py"
