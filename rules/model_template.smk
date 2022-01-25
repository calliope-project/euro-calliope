""" All scripts to produce Calliope model configuration YAML files """
import glob
from pathlib import Path

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"
template_dir = f"{root_dir}templates/"
locations_template_dir = f"{template_dir}techs_and_locations/"

NON_MODEL_CONFIG_FILES = ["environment.yaml", "README.md"]

MISSING_LOCATION_SPECIFIC_FILES = [
    "build/model/techs_and_locations/continental/transmission/transmission-electricity-entsoe.yaml",
    "build/model/techs_and_locations/regional/transmission/transmission-electricity-entsoe.yaml"
]


rule locations_template:
    message: "Generate locations configuration file for {wildcards.resolution} resolution from template."
    input:
        script = script_dir + "template_locations.py",
        template = locations_template_dir + "locations.yaml",
        shapes = rules.units.output[0]
    output:
        yaml = "build/model/techs_and_locations/{resolution}/locations.yaml",
        csv = "build/model/techs_and_locations/{resolution}/locations.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/template_locations.py"


rule simple_techs_and_locations_template:
    message: "Create {wildcards.resolution} tech definition file `{wildcards.template}` from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = locations_template_dir + "{template}",
        locations = rules.locations_template.output.csv
    params:
        scaling_factors = config["scaling-factors"],
    wildcard_constraints:
        template = "demand/demand-electricity.yaml|supply/supply-load-shedding.yaml|storage/storage-electricity.yaml"
    conda: "../envs/default.yaml"
    output: "build/model/techs_and_locations/{resolution}/{template}"
    script: "../scripts/template_techs_and_locations.py"


rule bio_techs_and_locations_template:
    message: "Create biofuel tech definition file from template."
    input:
        script = script_dir + "biofuels/template_bio.py",
        template = locations_template_dir + "supply/supply-biofuel.yaml",
        biofuel_cost = "build/data/regional/biofuel/{scenario}/costs-eur-per-mwh.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        ),
        locations = "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
    params:
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"],
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/model/techs_and_locations/{resolution}/supply/supply-biofuel.yaml"
    script: "../scripts/biofuels/template_bio.py"


rule hydro_supply_techs_at_locations_template:
    message: "Allocate hydro supply techs to {wildcards.resolution} locations from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = locations_template_dir + "supply/supply-hydro.yaml",
        locations = "build/data/{resolution}/hydro-capacities-mw.csv"
    params:
        capacity_factors = config["capacity-factors"]["average"],
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/model/techs_and_locations/{resolution}/supply/supply-hydro.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule hydro_storage_techs_at_locations_template:
    message: "Allocate hydro storage techs to {wildcards.resolution} locations from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = locations_template_dir + "storage/storage-hydro-electricity.yaml",
        locations = "build/data/{resolution}/hydro-capacities-mw.csv"
    params:
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/model/techs_and_locations/{resolution}/storage/storage-hydro-electricity.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule wind_solar_techs_at_locations_template:
    message: "Create {wildcards.resolution} wind & solar tech definition file from template."
    input:
        script = script_dir + "template_techs_and_locations.py",
        template = locations_template_dir + "supply/supply-wind-and-solar.yaml",
        locations = rules.area_to_capacity_limits.output[0]
    params:
        capacity_factors = config["capacity-factors"]["average"],
        scaling_factors = config["scaling-factors"],
        max_power_density = config["parameters"]["maximum-installable-power-density"]
    conda: "../envs/default.yaml"
    output: "build/model/techs_and_locations/{resolution}/supply/supply-wind-and-solar.yaml"
    script: "../scripts/template_techs_and_locations.py"


rule transmission_entsoe_tyndp_template:
    message: "Create YAML file of national-scale links with ENTSO-E TYNDP net-transfer capacities"
    input:
        script = script_dir + "transmission/template_transmission_entsoe_tyndp.py",
        template = locations_template_dir + "transmission/transmission-electricity-entsoe.yaml",
        locations = rules.locations_template.output.csv,
        entsoe_tyndp = "build/data/national/TYNDP-2020-Scenario-Datafile.xlsx"
    params:
        scenario = config["parameters"]["entsoe-tyndp"]["scenario"],
        grid = config["parameters"]["entsoe-tyndp"]["grid"],
        ntc_limit = config["parameters"]["entsoe-tyndp"]["ntc_limit"],
        energy_cap_limit = config["parameters"]["entsoe-tyndp"]["energy_cap_limit"],
        year = config["parameters"]["entsoe-tyndp"]["projection-year"],
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/model/techs_and_locations/{resolution}/transmission/transmission-electricity-entsoe.yaml"
    wildcard_constraints: resolution = "national"
    conda: "../envs/default.yaml"
    script: "../scripts/transmission/template_transmission_entsoe_tyndp.py"


rule link_locations_with_transmission_techs_template:
    message: "Link {wildcards.resolution} direct neighbours and neighbours with sea connections with transmission techs from template."
    input:
        script = script_dir + "transmission/template_transmission.py",
        template = locations_template_dir + "transmission/transmission-electricity.yaml",
        units = rules.units.output[0]
    params:
        scaling_factors = config["scaling-factors"],
        sea_connections = lambda wildcards: config["sea-connections"][wildcards.resolution]
    output: "build/model/techs_and_locations/{resolution}/transmission/transmission-electricity.yaml"
    conda: "../envs/geo.yaml"
    script: "../scripts/transmission/template_transmission.py"


rule no_params_template:
    message: "Create files from templates where no parameterisation is required."
    input:
        template = template_dir + "{template}",
    output: "build/model/{template}"
    wildcard_constraints:
        template = "README.md|environment.yaml|interest-rate.yaml|scenarios.yaml"
    conda: "../envs/shell.yaml"
    shell: "cp {input.template} {output}"


def config_files(other_model_files, resolution, missing_location_specific_files):
    top_level_files = [
        i.replace("templates/", "build/model/")
        for i in glob.glob("templates/*.yaml")
        if Path(i).name not in other_model_files
    ]
    locations_files = [
        i.replace("templates/techs_and_locations/", f"build/model/techs_and_locations/{resolution}/")
        for i in glob.glob("templates/techs_and_locations/**/*.yaml", recursive=True)
    ]
    locations_files = set(locations_files).difference(missing_location_specific_files)
    return [*top_level_files, *locations_files]


rule model_template:
    message: "Generate top-level {wildcards.resolution} model configuration file from template"
    input:
        script = script_dir + "template_model.py",
        template = template_dir + "example-model.yaml",
        other_model_files = expand(
            "build/model/{template}", template=NON_MODEL_CONFIG_FILES
        ),
        config_files = lambda wildcards: config_files(
            NON_MODEL_CONFIG_FILES + ["example-model.yaml"],
            wildcards.resolution,
            MISSING_LOCATION_SPECIFIC_FILES
        )
    params:
        model_year = config["parameters"]["model-year"]
    conda: "../envs/default.yaml"
    output: "build/model/example-model-{resolution}.yaml"
    script: "../scripts/template_model.py"
