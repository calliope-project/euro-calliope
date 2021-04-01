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
localrules: raw_gadm_administrative_borders_zipped, raw_gadm_administrative_borders, raw_nuts_units_zipped
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
src_dir = f"{root_dir}src/"


rule raw_gadm_administrative_borders_zipped:
    message: "Download administrative borders for {wildcards.country_code} as zip."
    params: url = config["data-sources"]["gadm"]
    output: protected("data/automatic/raw-gadm/{country_code}.zip")
    shell: "curl -sLo {output} '{params.url}{wildcards.country_code}_gpkg.zip'"


rule raw_gadm_administrative_borders:
    message: "Unzip administrative borders of {wildcards.country_code} as zip."
    input: "data/automatic/raw-gadm/{country_code}.zip"
    output: temp("data/automatic/raw-gadm/gadm36_{country_code}.gpkg")
    shell: "unzip -o {input} -d data/automatic/raw-gadm"


rule administrative_borders_gadm:
    message: "Merge administrative borders of all countries up to layer {params.max_layer_depth}."
    input:
        src = src_dir + "shapes/gadm.py",
        countries = ["data/automatic/raw-gadm/gadm36_{}.gpkg".format(country_code)
                     for country_code in [pycountry.countries.lookup(country).alpha_3
                                          for country in config['scope']['countries']]
                    ]
    params:
        max_layer_depth = 2,
        crs = config["crs"],
        schema = SCHEMA_UNITS,
        x_min = config["scope"]["bounds"]["x_min"],
        x_max = config["scope"]["bounds"]["x_max"],
        y_min = config["scope"]["bounds"]["y_min"],
        y_max = config["scope"]["bounds"]["y_max"]
    output: "build/data/administrative-borders-gadm.gpkg"
    conda: "../envs/geo.yaml"
    script: "../src/shapes/gadm.py"


rule raw_nuts_units_zipped:
    message: "Download units as zip."
    params: url = config["data-sources"]["nuts"]
    output: protected("data/automatic/raw-nuts-units.zip")
    shell: "curl -sLo {output} '{params.url}'"


rule administrative_borders_nuts:
    message: "Normalise NUTS administrative borders."
    input:
        src = src_dir + "shapes/nuts.py",
        zipped = rules.raw_nuts_units_zipped.output[0]
    params:
        crs = config["crs"],
        schema = SCHEMA_UNITS,
        x_min = config["scope"]["bounds"]["x_min"],
        x_max = config["scope"]["bounds"]["x_max"],
        y_min = config["scope"]["bounds"]["y_min"],
        y_max = config["scope"]["bounds"]["y_max"],
        all_countries = config["scope"]["countries"]
    output: "build/data/administrative-borders-nuts.gpkg"
    shadow: "minimal"
    conda: "../envs/geo.yaml"
    script: "../src/shapes/nuts.py"


rule units:
    message: "Form units of resolution {wildcards.resolution} by remixing NUTS and GADM."
    input:
        src = src_dir + "shapes/units.py",
        nuts = rules.administrative_borders_nuts.output[0],
        gadm = rules.administrative_borders_gadm.output[0]
    params:
        all_countries = config["scope"]["countries"],
        layer_configs = config["shapes"]
    output:
        "build/data/{resolution}/units.geojson"
    conda: "../envs/geo.yaml"
    script: "../src/shapes/units.py"


rule units_without_shape:
    message: "Dataset of units on resolution {wildcards.resolution} without geo information."
    input:
        src = src_dir + "shapes/nogeo.py",
        units = rules.units.output[0]
    output: "build/data/{resolution}/units.csv"
    conda: "../envs/geo.yaml"
    script: "../src/shapes/nogeo.py"


rule eez:
    message: "Clip exclusive economic zones to study area."
    input: config["data-sources"]["eez"]
    output: "build/data/eez.geojson"
    params:
        bounds="{x_min},{y_min},{x_max},{y_max}".format(**config["scope"]["bounds"]),
        countries=",".join(["'{}'".format(country) for country in config["scope"]["countries"]])
    conda: "../envs/geo.yaml"
    shell:
        """
        fio cat --bbox {params.bounds} {input}\
        | fio filter "f.properties.Territory1 in [{params.countries}]"\
        | fio collect > {output}
        """
