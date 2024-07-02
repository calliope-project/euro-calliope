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
include: "./rules/jrc-idees.smk"
include: "./rules/wind-and-solar.smk"
include: "./rules/biofuels.smk"
include: "./rules/hydro.smk"
include: "./rules/transmission.smk"
include: "./rules/demand.smk"
include: "./rules/nuclear.smk"
include: "./rules/transport.smk"
include: "./rules/sync.smk"
include: "./rules/heat.smk"
include: "./rules/modules.smk"
min_version("8.10")
localrules: all, clean
wildcard_constraints:
    resolution = "continental|national|regional|ehighways"

ruleorder: techs_and_locations_template > dummy_tech_locations_template

ALL_CF_TECHNOLOGIES = [
    "wind-onshore", "wind-offshore", "open-field-pv",
    "rooftop-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat", "hydro-run-of-river",
    "hydro-reservoir"
]


def ensure_lib_folder_is_linked():
    if not (hasattr(workflow, "deployment_settings") and not
            hasattr(workflow.deployment_settings, "conda_prefix")):
        return
    link = Path(workflow.deployment_settings.conda_prefix) / "lib"
    if not link.exists():
        # Link either does not exist or is an invalid symlink
        print("Creating link from conda env dir to eurocalliopelib.")
        if link.is_symlink():  # Deal with existing but invalid symlink
            shell(f"rm {link}")
        makedirs(workflow.deployment_settings.conda_prefix)
        shell(f"ln -s {workflow.basedir}/lib {link}")


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
    localrule: True
    default_target: True
    input:
        "build/logs/continental/test.success",
        "build/logs/national/test.success",
        expand(
            "build/models/{resolution}/{file}",
            resolution=["continental", "national", "regional", "ehighways"],
            file=["example-model.yaml", "build-metadata.yaml", "summary-of-potentials.nc", "summary-of-potentials.csv"]
        )

rule all_tests:
    message: "Generate euro-calliope pre-built models and run all tests."
    input:
        expand(
            "build/logs/{resolution}/test.success",
            resolution=["continental", "national", "regional", "ehighways"],
        ),
        expand(
            "build/models/{resolution}/{file}",
            resolution=["continental", "national", "regional", "ehighways"],
            file=["example-model.yaml", "build-metadata.yaml", "summary-of-potentials.nc", "summary-of-potentials.csv"]
        )


rule techs_and_locations_template:
    message: "Create {wildcards.resolution} definition file for {wildcards.tech_and_group}."
    input:
        template = techs_template_dir + "{tech_and_group}.yaml.jinja",
        locations = "build/data/{resolution}/{tech_and_group}.csv"
    params:
        scaling_factors = config["scaling-factors"],
        capacity_factors = config["capacity-factors"]["average"],
        max_power_densities = config["parameters"]["maximum-installable-power-density"],
        heat_pump_shares = config["parameters"]["heat-pump"]["heat-pump-shares"],
    wildcard_constraints:
        # Exclude all outputs that have their own `techs_and_locations_template` implementation
        tech_and_group = "(?!transmission\/|supply\/biofuel).*"
    conda: "envs/default.yaml"
    output: "build/models/{resolution}/techs/{tech_and_group}.yaml"
    script: "scripts/template_techs.py"

use rule techs_and_locations_template as dummy_tech_locations_template with:
    # For all cases where we don't have any location-specific data that we want to supply to the template
    input:
        template = techs_template_dir + "{tech_and_group}.yaml.jinja",
        locations = rules.locations_template.output.csv

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
        template = model_template_dir + "example-model.yaml.jinja",
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
                "techs/demand/heat.yaml",
                "techs/storage/electricity.yaml",
                "techs/storage/hydro.yaml",
                "techs/storage/heat.yaml",
                "techs/supply/biofuel.yaml",
                "techs/supply/hydro.yaml",
                "techs/supply/load-shedding.yaml",
                "techs/supply/open-field-solar-and-wind-onshore.yaml",
                "techs/supply/rooftop-solar.yaml",
                "techs/supply/wind-offshore.yaml",
                "techs/supply/nuclear.yaml",
                "techs/supply/heat-from-electricity.yaml",
                "techs/supply/historic-electrified-heat.yaml",
            ]
        ),
        heat_supply = (
            "build/models/{resolution}/timeseries/supply/heat-pump-cop.csv",
            "build/models/{resolution}/timeseries/supply/historic-electrified-heat.csv",
        ),
        capacityfactor_timeseries_data = expand(
            "build/models/{{resolution}}/timeseries/supply/capacityfactors-{technology}.csv",
            technology=ALL_CF_TECHNOLOGIES
        ),
        demand_timeseries_data = (
            "build/models/{resolution}/timeseries/demand/electricity.csv",
            "build/models/{resolution}/timeseries/demand/uncontrolled-electrified-road-transport.csv",
            "build/models/{resolution}/timeseries/demand/uncontrolled-road-transport-historic-electrification.csv",
            "build/models/{resolution}/timeseries/demand/heat.csv",
            "build/models/{resolution}/timeseries/demand/electrified-heat.csv",
            "build/models/{resolution}/timeseries/demand/demand-shape-min-ev.csv",
            "build/models/{resolution}/timeseries/demand/demand-shape-max-ev.csv",
            "build/models/{resolution}/timeseries/demand/demand-shape-equals-ev.csv",
            "build/models/{resolution}/timeseries/demand/plugin-profiles-ev.csv",
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
    localrule: True
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
        ),
        electrified_heat_demand = "build/models/{resolution}/timeseries/demand/electrified-heat.csv",
        heat_demand = "build/models/{resolution}/timeseries/demand/heat.csv",
        historic_electrified_heat = "build/models/{resolution}/timeseries/supply/historic-electrified-heat.csv",
        cop = "build/models/{resolution}/timeseries/supply/heat-pump-cop.csv"
    params:
        config = config,
        test_args = []  # add e.g. "--pdb" to enter ipdb on test failure
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
