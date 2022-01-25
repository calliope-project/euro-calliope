""" All scripts to produce Calliope model configuration YAML files """
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"
template_dir = f"{root_dir}templates/"
locations_template_dir = f"{template_dir}techs_and_locations/"



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


