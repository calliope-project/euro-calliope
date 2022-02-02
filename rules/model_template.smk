""" All scripts to produce Calliope model configuration YAML files """
import glob
from pathlib import Path

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"
template_dir = f"{root_dir}templates/"
model_template_dir = f"{template_dir}model/"
techs_template_dir = f"{model_template_dir}techs/"

NON_MODEL_CONFIG_FILES = ["environment.yaml", "README.md"]

MISSING_LOCATION_SPECIFIC_FILES = {
    "continental": ["techs/transmission/transmission-electricity-entsoe.yaml"],
    "regional": ["techs/transmission/transmission-electricity-entsoe.yaml"]
}


rule locations_template:
    message: "Generate locations configuration file for {wildcards.resolution} resolution from template."
    input:
        script = script_dir + "template_locations.py",
        template = model_template_dir + "locations.yaml",
        shapes = rules.units.output[0]
    output:
        yaml = "build/models/{resolution}/locations.yaml",
        csv = "build/data/{resolution}/locations.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/template_locations.py"


rule simple_techs_and_locations_template:
    message: "Create {wildcards.resolution} tech definition file `{wildcards.template}` from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = techs_template_dir + "{template}",
        locations = rules.locations_template.output.csv
    params:
        scaling_factors = config["scaling-factors"],
    wildcard_constraints:
        template = "demand/demand-electricity.yaml|supply/supply-load-shedding.yaml|storage/storage-electricity.yaml"
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/{template}"
    script: "../scripts/template_techs_and_locations.py"


rule bio_techs_and_locations_template:
    message: "Create biofuel tech definition file from template."
    input:
        script = script_dir + "biofuels/template_bio.py",
        template = techs_template_dir + "supply/supply-biofuel.yaml",
        biofuel_cost = "build/data/regional/biofuel/{scenario}/costs-eur-per-mwh.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        ),
        locations = "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
    params:
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"],
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/supply/supply-biofuel.yaml"
    script: "../scripts/biofuels/template_bio.py"


rule hydro_supply_techs_at_locations_template:
    message: "Allocate hydro supply techs to {wildcards.resolution} locations from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = techs_template_dir + "supply/supply-hydro.yaml",
        locations = "build/data/{resolution}/hydro-capacities-mw.csv"
    params:
        capacity_factors = config["capacity-factors"]["average"],
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/supply/supply-hydro.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule hydro_storage_techs_at_locations_template:
    message: "Allocate hydro storage techs to {wildcards.resolution} locations from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = techs_template_dir + "storage/storage-hydro-electricity.yaml",
        locations = "build/data/{resolution}/hydro-capacities-mw.csv"
    params:
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/storage/storage-hydro-electricity.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule wind_solar_techs_at_locations_template:
    message: "Create {wildcards.resolution} wind & solar tech definition file from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = techs_template_dir + "supply/supply-wind-and-solar.yaml",
        locations = rules.area_to_capacity_limits.output[0]
    params:
        capacity_factors = config["capacity-factors"]["average"],
        scaling_factors = config["scaling-factors"],
        max_power_density = config["parameters"]["maximum-installable-power-density"]
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/supply/supply-wind-and-solar.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule transmission_entsoe_tyndp_template:
    message: "Create YAML file of national-scale links with ENTSO-E TYNDP net-transfer capacities"
    input:
        script = script_dir + "transmission/template_transmission_entsoe_tyndp.py",
        template = techs_template_dir + "transmission/transmission-electricity-entsoe.yaml",
        locations = rules.locations_template.output.csv,
        entsoe_tyndp = "build/data/national/TYNDP-2020-Scenario-Datafile.xlsx"
    params:
        scenario = config["parameters"]["entsoe-tyndp"]["scenario"],
        grid = config["parameters"]["entsoe-tyndp"]["grid"],
        ntc_limit = config["parameters"]["entsoe-tyndp"]["ntc_limit"],
        energy_cap_limit = config["parameters"]["entsoe-tyndp"]["energy_cap_limit"],
        year = config["parameters"]["entsoe-tyndp"]["projection-year"],
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/techs/transmission/transmission-electricity-entsoe.yaml"
    wildcard_constraints: resolution = "national"
    conda: "../envs/default.yaml"
    script: "../scripts/transmission/template_transmission_entsoe_tyndp.py"


rule link_locations_with_transmission_techs_template:
    message: "Link {wildcards.resolution} direct neighbours and neighbours with sea connections with transmission techs from template."
    input:
        script = script_dir + "transmission/template_transmission.py",
        template = techs_template_dir + "transmission/transmission-electricity.yaml",
        units = rules.units.output[0]
    params:
        scaling_factors = config["scaling-factors"],
        sea_connections = lambda wildcards: config["sea-connections"][wildcards.resolution]
    output: "build/models/{resolution}/techs/transmission/transmission-electricity.yaml"
    conda: "../envs/geo.yaml"
    script: "../scripts/transmission/template_transmission.py"


rule no_params_model_template:
    message: "Create {wildcards.resolution} configuration files from templates where no parameterisation is required."
    input:
        template = model_template_dir + "{template}",
    output: "build/models/{resolution}/{template}"
    wildcard_constraints:
        template = "interest-rate.yaml|scenarios.yaml"
    conda: "../envs/shell.yaml"
    shell: "cp {input.template} {output}"


rule no_params_template:
    message: "Create non-model files from templates where no parameterisation is required."
    input:
        template = template_dir + "{template}",
    output: "build/models/{template}"
    wildcard_constraints:
        template = "|".join(NON_MODEL_CONFIG_FILES)
    conda: "../envs/shell.yaml"
    shell: "cp {input.template} {output}"


def get_all_config_files(resolution):
    all_template_files = glob.glob("templates/model/**/*.yaml", recursive=True)
    all_files = set([i.replace("templates/model/", "") for i in all_template_files])
    all_files_strip_missing = all_files.difference(
        MISSING_LOCATION_SPECIFIC_FILES.get(resolution, []) + ["example-model.yaml"]
    )
    return list(all_files_strip_missing)


rule model_template:
    message: "Generate top-level {wildcards.resolution} model configuration file from template"
    input:
        script = script_dir + "template_model.py",
        template = model_template_dir + "example-model.yaml",
        non_model_files = expand(
            "build/models/{template}", template=NON_MODEL_CONFIG_FILES
        ),
        config_files = lambda wildcards: expand(
            f"build/models/{wildcards.resolution}/{{filepath}}",
            filepath=get_all_config_files(wildcards.resolution)
        )
    params:
        model_year = config["parameters"]["model-year"]
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/example-model.yaml"
    script: "../scripts/template_model.py"
