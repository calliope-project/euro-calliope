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
