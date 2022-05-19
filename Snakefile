import glob
from pathlib import Path

from snakemake.utils import validate

configfile: "config/default.yaml"
validate(config, "config/schema.yaml")

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
__version__ = open(f"{root_dir}VERSION").readlines()[0].strip()
script_dir = f"{root_dir}scripts/"
test_dir = f"{root_dir}tests/"
model_test_dir = f"{test_dir}model"
template_dir = f"{root_dir}templates/"
model_template_dir = f"{template_dir}models/"
techs_template_dir = f"{model_template_dir}techs/"

include: "./rules/shapes.smk"
include: "./rules/wind-and-solar.smk"
include: "./rules/biofuels.smk"
include: "./rules/hydro.smk"
include: "./rules/transmission.smk"
include: "./rules/demand.smk"
include: "./rules/nuclear.smk"
include: "./rules/sync.smk"
localrules: all, clean
wildcard_constraints:
        resolution = "continental|national|regional"

ruleorder: area_to_capacity_limits > hydro_capacities > biofuels > nuclear_regional_capacity > dummy_tech_locations_template
ruleorder: bio_techs_and_locations_template > techs_and_locations_template

ALL_CF_TECHNOLOGIES = [
    "wind-onshore", "wind-offshore", "open-field-pv",
    "rooftop-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat", "hydro-run-of-river",
    "hydro-reservoir"
]
ALL_DEMAND_CARRIERS = ["electricity"]

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
        "build/logs/continental/test.success",
        "build/logs/national/test.success",
        "build/models/continental/example-model.yaml",
        "build/models/national/example-model.yaml",
        "build/models/regional/example-model.yaml",
        "build/models/build-metadata.yaml"


rule all_tests:
    message: "Generate euro-calliope pre-built models and run all tests."
    input:
        "build/models/continental/example-model.yaml",
        "build/models/national/example-model.yaml",
        "build/models/regional/example-model.yaml",
        "build/logs/continental/test.success",
        "build/logs/national/test.success",
        "build/logs/regional/test.success",
        "build/models/build-metadata.yaml"


rule dummy_tech_locations_template:  # needed to provide `techs_and_locations_template` with a locational CSV linked to each technology that has no location-specific data to define.
    message: "Create empty {wildcards.resolution} location-specific data file for the {wildcards.tech_group} tech `{wildcards.tech}`."
    input: rules.locations_template.output.csv
    output: "build/data/{resolution}/{tech_group}/{tech}.csv"
    conda: "envs/shell.yaml"
    shell: "cp {input} {output}"


rule techs_and_locations_template:
    message: "Create {wildcards.resolution} definition file for the {wildcards.tech_group} tech `{wildcards.tech}`."
    input:
        script = script_dir + "template_techs.py",
        template = techs_template_dir + "{tech_group}/{tech}.yaml",
        locations = "build/data/{resolution}/{tech_group}/{tech}.csv"
    params:
        scaling_factors = config["scaling-factors"],
        capacity_factors = config["capacity-factors"]["average"],
        max_power_densities = config["parameters"]["maximum-installable-power-density"]
    wildcard_constraints:
        tech_group = "(?!transmission).*"  # i.e. all but transmission
    conda: "envs/default.yaml"
    output: "build/models/{resolution}/techs/{tech_group}/{tech}.yaml"
    script: "scripts/template_techs.py"


rule no_params_model_template:
    message: "Create {wildcards.resolution} configuration files from templates where no parameterisation is required."
    input:
        template = model_template_dir + "{template}",
    output: "build/models/{resolution}/{template}"
    wildcard_constraints:
        template = "interest-rate.yaml|scenarios.yaml"
    conda: "envs/shell.yaml"
    shell: "cp {input.template} {output}"


rule no_params_template:
    message: "Create non-model files from templates where no parameterisation is required."
    input:
        template = template_dir + "{template}",
    output: "build/models/{template}"
    wildcard_constraints:
        template = "[^/]*"
    conda: "envs/shell.yaml"
    shell: "cp {input.template} {output}"


rule model_template:
    message: "Generate top-level {wildcards.resolution} model configuration file from template"
    input:
        script = script_dir + "template_model.py",
        template = model_template_dir + "example-model.yaml",
        non_model_files = expand(
            "build/models/{template}", template=["environment.yaml", "README.md"]
        ),
        input_files = expand(
            "build/models/{{resolution}}/{input_file}",
            input_file=[
                "interest-rate.yaml",
                "locations.yaml",
                "scenarios.yaml",
                "techs/demand/electricity.yaml",
                "techs/storage/electricity.yaml",
                "techs/storage/hydro.yaml",
                "techs/supply/biofuel.yaml",
                "techs/supply/hydro.yaml",
                "techs/supply/load-shedding.yaml",
                "techs/supply/open-field-solar-and-wind-onshore.yaml",
                "techs/supply/rooftop-solar.yaml",
                "techs/supply/wind-offshore.yaml",
                "techs/supply/nuclear.yaml",
            ]
        ),
        capacityfactor_timeseries_data = expand(
            "build/models/{{resolution}}/timeseries/supply/capacityfactors-{technology}.csv",
            technology=ALL_CF_TECHNOLOGIES
        ),
        demand_timeseries_data = expand(
            "build/models/{{resolution}}/timeseries/demand/{energy_carrier}.csv",
            energy_carrier=ALL_DEMAND_CARRIERS
        ),
        optional_input_files = lambda wildcards: expand(
            f"build/models/{wildcards.resolution}/{{input_file}}",
            input_file=[
                "techs/transmission/electricity-linked-neighbours.yaml",
            ] + ["techs/transmission/electricity-entsoe.yaml" for i in [None] if wildcards.resolution == "national"]
        )

    params:
        year = config["scope"]["temporal"]["first-year"]
    conda: "envs/default.yaml"
    output: "build/models/{resolution}/example-model.yaml"
    script: "scripts/template_model.py"


rule build_metadata:
    message: "Generate build metadata."
    input:
        script_dir + "metadata.py",
        "build/models/continental/example-model.yaml",
        "build/models/national/example-model.yaml",
        "build/models/regional/example-model.yaml",
    params:
        config = config,
        version = __version__
    output: "build/models/build-metadata.yaml"
    conda: "envs/default.yaml"
    script: "scripts/metadata.py"


rule clean: # removes all generated results
    shell:
        """
        rm -r build/
        echo "Data downloaded to data/automatic/ has not been cleaned."
        """


rule test:
    message: "Run tests"
    input:
        test_dir = model_test_dir,
        tests = map(str, Path(model_test_dir).glob("**/test_*.py")),
        example_model = "build/models/{resolution}/example-model.yaml",
        capacity_factor_timeseries = expand(
            "build/models/{{resolution}}/timeseries/supply/capacityfactors-{technology}.csv",
            technology=ALL_CF_TECHNOLOGIES
        )
    params:
        config = config
    log: "build/logs/{resolution}/test-report.html"
    output: "build/logs/{resolution}/test.success"
    conda: "./envs/test.yaml"
    script: "./tests/model/test_runner.py"
