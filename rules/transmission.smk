"""Rules to generate electricity transmission lines between regions."""

localrules: download_entsoe_tyndp_zip, entsoe_tyndp_xlsx


rule download_entsoe_tyndp_zip:
    message: "Download ENTSO-E ten-year network development plan (TYNDP) 2020 scenario dataset"
    params: url = config["data-sources"]["entsoe-tyndp"]
    output: protected("data/automatic/raw-entsoe-tyndp.xlsx.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule entsoe_tyndp_xlsx:
    message: "Unzip ENTSO-E TYNDP 2020 scenario dataset."
    input: rules.download_entsoe_tyndp_zip.output[0]
    shadow: "minimal"
    output: "build/data/national/TYNDP-2020-Scenario-Datafile.xlsx",
    conda: "../envs/shell.yaml"
    shell: "unzip -o {input} 'TYNDP-2020-Scenario-Datafile.xlsx' -d build/data/national"


rule transmission_entsoe_tyndp_template:
    message: "Create YAML file of national-scale links with ENTSO-E TYNDP net-transfer capacities"
    input:
        script = script_dir + "transmission/template_transmission_entsoe_tyndp.py",
        template = techs_template_dir + "transmission/electricity-transmission.yaml",
        locations = rules.locations_template.output.csv,
        entsoe_tyndp = rules.entsoe_tyndp_xlsx.output[0]
    params:
        scenario = config["parameters"]["entsoe-tyndp"]["scenario"],
        grid = config["parameters"]["entsoe-tyndp"]["grid"],
        ntc_limit = config["parameters"]["entsoe-tyndp"]["ntc_limit"],
        energy_cap_limit = config["parameters"]["entsoe-tyndp"]["energy_cap_limit"],
        year = config["parameters"]["entsoe-tyndp"]["projection-year"],
        scaling_factors = config["scaling-factors"]
    output: "build/models/{resolution}/techs/transmission/electricity-entsoe.yaml"
    wildcard_constraints: resolution = "national"
    conda: "../envs/default.yaml"
    script: "../scripts/transmission/template_transmission_entsoe_tyndp.py"


rule link_locations_with_transmission_techs_template:
    message: "Link {wildcards.resolution} direct neighbours and neighbours with sea connections with transmission techs from template."
    input:
        script = script_dir + "transmission/template_transmission.py",
        template = techs_template_dir + "transmission/electricity-transmission.yaml",
        units = rules.units.output[0]
    params:
        scaling_factors = config["scaling-factors"],
        sea_connections = lambda wildcards: config["sea-connections"][wildcards.resolution]
    output: "build/models/{resolution}/techs/transmission/electricity-linked-neighbours.yaml"
    conda: "../envs/geo.yaml"
    script: "../scripts/transmission/template_transmission.py"
