"""Snakemake rules to create shapes of administrative regions and exclusive economic zones."""
import pycountry

SCHEMA_UNITS = {
    "properties": {
        "country_code": "str", # id of the country to which the unit belongs
        "id": "str", # a unique id of this unit
        "name": "str", # the name of the unit, not necessarily unqique
        "type": "str", # the type of the unit
        "proper": "bool" # flag indicating proper administrative unit (not the case for water bodies e.g.)
    },
    "geometry": "MultiPolygon"
}

configfile: "config/default.yaml"
localrules: download_raw_gadm_administrative_borders, raw_gadm_administrative_borders, download_raw_nuts_units
localrules: download_eez
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule download_raw_gadm_administrative_borders:
    message: "Download administrative borders for {wildcards.country_code} as zip."
    params: url = lambda wildcards: config["data-sources"]["gadm"].format(country_code=wildcards.country_code)
    output: protected("data/automatic/raw-gadm/{country_code}.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule raw_gadm_administrative_borders:
    message: "Unzip administrative borders of {wildcards.country_code} as zip."
    input: rules.download_raw_gadm_administrative_borders.output[0]
    output: temp("build/data/raw-gadm/gadm36_{country_code}.gpkg")
    conda: "../envs/shell.yaml"
    shell: "unzip -o {input} -d build/data/raw-gadm"


rule administrative_borders_gadm:
    message: "Merge administrative borders of all countries up to layer {params.max_layer_depth}."
    input:
        script = script_dir + "shapes/gadm.py",
        countries = [f"build/data/raw-gadm/gadm36_{country_code}.gpkg"
                     for country_code in [pycountry.countries.lookup(country).alpha_3
                                          for country in config["scope"]["spatial"]["countries"]]
                    ]
    params:
        max_layer_depth = 2,
        crs = config["crs"],
        schema = SCHEMA_UNITS,
        x_min = config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = config["scope"]["spatial"]["bounds"]["y_max"]
    output: "build/data/administrative-borders-gadm.gpkg"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/gadm.py"


rule download_raw_nuts_units:
    message: "Download units as zip."
    params: url = config["data-sources"]["nuts"]
    output: protected("data/automatic/raw-nuts-units.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule administrative_borders_nuts:
    message: "Normalise NUTS administrative borders."
    input:
        script = script_dir + "shapes/nuts.py",
        zipped = rules.download_raw_nuts_units.output[0]
    params:
        crs = config["crs"],
        schema = SCHEMA_UNITS,
        x_min = config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = config["scope"]["spatial"]["bounds"]["y_max"],
        all_countries = config["scope"]["spatial"]["countries"]
    output: "build/data/administrative-borders-nuts.gpkg"
    shadow: "minimal"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/nuts.py"


rule units:
    message: "Form units of resolution {wildcards.resolution} by remixing NUTS and GADM."
    input:
        script = script_dir + "shapes/units.py",
        nuts = rules.administrative_borders_nuts.output[0],
        gadm = rules.administrative_borders_gadm.output[0]
    params:
        all_countries = config["scope"]["spatial"]["countries"],
        layer_configs = config["shapes"]
    output:
        "build/data/{resolution}/units.geojson"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/units.py"


rule units_without_shape:
    message: "Dataset of units on resolution {wildcards.resolution} without geo information."
    input:
        script = script_dir + "shapes/nogeo.py",
        units = rules.units.output[0]
    output: "build/data/{resolution}/units.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/nogeo.py"


rule download_eez:
    message: "Download Exclusive Economic Zones as zip"
    output: protected("data/automatic/eez.zip")
    params: url = config["data-sources"]["eez"]
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule eez:
    message: "Clip exclusive economic zones to study area."
    input: rules.download_eez.output[0]
    output: "build/data/eez.geojson"
    params:
        bounds="{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["spatial"]["bounds"]),
        countries=",".join(["'{}'".format(country) for country in config["scope"]["spatial"]["countries"]]),
    conda: "../envs/geo.yaml"
    shadow: "minimal"
    shell:
        """
        fio cat --bbox {params.bounds} "zip://{input}"\
        | fio filter "f.properties.territory1 in [{params.countries}]"\
        | fio collect > {output}
        """
