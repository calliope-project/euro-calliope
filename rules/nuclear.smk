localrules: jrc_power_plant_database_zipped

rule jrc_power_plant_database_zipped:
    message: "Download and unzip jrc power plant database."
    params: url = config["data-sources"]["jrc-ppdb"]
    output:  "data/automatic/jrc_power_plant_database.zip"
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output[0]} '{params.url}'"


rule jrc_power_plant_database:
    message: "Unzip JRC power plant units from database."
    input: rules.jrc_power_plant_database_zipped.output[0]
    output:
        "data/automatic/JRC_OPEN_UNITS.csv"
    conda: "../envs/shell.yaml"
    shell: "unzip -o {input[0]} JRC_OPEN_UNITS.csv -d data/automatic/"


rule nuclear_regional_capacity:
    message: "Use current geolocations of nuclear capacity in Europe to assign nuclear capacity to {wildcards.resolution} regions."
    input:
        script = script_dir + "nuclear/regional_capacity.py",
        power_plant_database = rules.jrc_power_plant_database.output[0],
        units = rules.units.output[0]
    params:
        nuclear_capacity_scenario = config["parameters"]["nuclear-capacity-scenario"],
        countries = config["scope"]["spatial"]["countries"]
    conda: "../envs/geo.yaml"
    output: "build/data/{resolution}/supply/nuclear.csv"
    script: "../scripts/nuclear/regional_capacity.py"
