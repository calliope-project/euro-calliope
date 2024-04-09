import glob
from pathlib import Path

from snakemake.utils import validate, min_version, makedirs

configfile: "config/default.yaml"
validate(config, "config/schema.yaml")

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
__version__ = open(f"{root_dir}VERSION").readlines()[0].strip()
test_dir = f"{root_dir}tests/"
model_test_dir = f"{test_dir}model"
template_dir = f"{root_dir}templates/"
model_template_dir = f"{template_dir}models/"
techs_template_dir = f"{model_template_dir}techs/"

include: "./rules/shapes.smk"
include: "./rules/data.smk"
include: "./rules/wind-and-solar.smk"
include: "./rules/biofuels.smk"
include: "./rules/hydro.smk"
include: "./rules/transmission.smk"
include: "./rules/demand.smk"
include: "./rules/nuclear.smk"
include: "./rules/transport.smk"
include: "./rules/sync.smk"
include: "./rules/heat.smk"
min_version("7.8")
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


def ensure_lib_folder_is_linked():
    if not workflow.deployment_settings.conda_prefix:
        return
    link = Path(workflow.deployment_settings.conda_prefix) / "lib"
    if not link.exists():
        print("Creating link from conda env dir to eurocalliopelib.")
        makedirs(workflow.deployment_settings.conda_prefix)
        shell(f"ln -s {workflow.basedir}/lib {workflow.deployment_settings.conda_prefix}/lib")


ensure_lib_folder_is_linked()

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
        "build/models/continental/build-metadata.yaml",
        "build/models/national/build-metadata.yaml",
        "build/models/regional/build-metadata.yaml",
        "build/models/regional/summary-of-potentials.nc",
        "build/models/regional/summary-of-potentials.csv",
        "build/models/national/summary-of-potentials.nc",
        "build/models/national/summary-of-potentials.csv",
        "build/models/continental/summary-of-potentials.nc",
        "build/models/continental/summary-of-potentials.csv"


rule all_tests:
    message: "Generate euro-calliope pre-built models and run all tests."
    input:
        "build/models/continental/example-model.yaml",
        "build/models/national/example-model.yaml",
        "build/models/regional/example-model.yaml",
        "build/logs/continental/test.success",
        "build/logs/national/test.success",
        "build/logs/regional/test.success",
        "build/models/build-metadata.yaml",
        "build/models/regional/summary-of-potentials.nc",
        "build/models/regional/summary-of-potentials.csv",
        "build/models/national/summary-of-potentials.nc",
        "build/models/national/summary-of-potentials.csv",
        "build/models/continental/summary-of-potentials.nc",
        "build/models/continental/summary-of-potentials.csv"


rule dummy_tech_locations_template:  # needed to provide `techs_and_locations_template` with a locational CSV linked to each technology that has no location-specific data to define.
    message: "Create empty {wildcards.resolution} location-specific data file for the {wildcards.tech_group} tech `{wildcards.tech}`."
    input: rules.locations_template.output.csv
    output: "build/data/{resolution}/{tech_group}/{tech}.csv"
    conda: "envs/shell.yaml"
    shell: "cp {input} {output}"


rule techs_and_locations_template:
    message: "Create {wildcards.resolution} definition file for the {wildcards.tech_group} tech `{wildcards.tech}`."
    input:
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
        template = "interest-rate.yaml"
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
        template = model_template_dir + "example-model.yaml",
        non_model_files = expand(
            "build/models/{template}", template=["environment.yaml", "README.md"]
        ),
        input_files = expand(
            "build/models/{{resolution}}/{input_file}",
            input_file=[
                "interest-rate.yaml",
                "locations.yaml",
                "techs/demand/electricity.yaml",
                "techs/demand/electrified-transport.yaml",
                "techs/demand/electrified-heat.yaml",
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
        demand_timeseries_data = (
            "build/models/{resolution}/timeseries/demand/electricity.csv",
            "build/models/{resolution}/timeseries/demand/electrified-road-transport.csv",
            "build/models/{resolution}/timeseries/demand/road-transport-historic-electrification.csv",
            "build/models/{resolution}/timeseries/demand/electrified-heat-demand.csv",
            "build/models/{resolution}/timeseries/demand/heat-demand-historic-electrification.csv",
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
        "build/models/{resolution}/example-model.yaml",
    params:
        config = config,
        version = __version__
    output: "build/models/{resolution}/build-metadata.yaml"
    conda: "envs/default.yaml"
    script: "scripts/metadata.py"


rule dag_dot:
    output: temp("build/dag.dot")
    shell:
        "snakemake --rulegraph > {output}"


rule dag:
    message: "Plot dependency graph of the workflow."
    input: rules.dag_dot.output[0]
    # Output is deliberatly omitted so rule is executed each time.
    conda: "envs/dag.yaml"
    shell:
        "dot -Tpdf {input} -o build/dag.pdf"


rule clean:  # removes all generated results
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
    resources:
        runtime = 240
    script: "./tests/model/test_runner.py"


rule summarise_potentials:
    message: "Generates netcdf and csv file with potentials for each technology."
    input:
        path_to_model = "build/models/{resolution}/example-model.yaml"
    output:
        netcdf = "build/models/{resolution}/summary-of-potentials.nc",
        csv = "build/models/{resolution}/summary-of-potentials.csv"
    params:
        scaling_factors = config["scaling-factors"]
    conda:
        "./envs/test.yaml"
    script:
        "./scripts/summarise_potentials.py"
