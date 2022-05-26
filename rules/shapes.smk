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

localrules: download_raw_gadm_administrative_borders, raw_gadm_administrative_borders, download_raw_nuts_units_zip, download_raw_nuts_units_geojson
localrules: download_eez, download_industry_sites_zip


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


rule download_raw_nuts_units_zip:
    message: "Download units as zip."
    params: url = config["data-sources"]["nuts"]
    output: protected("data/automatic/raw-nuts-units.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule download_raw_nuts_units_geojson:
    message: "Download NUTS{wildcards.nuts_year} units as geojson."
    params: url = lambda wildcards: config["data-sources"]["nuts-geojson"].format(year=wildcards.nuts_year)
    output: protected("data/automatic/raw-nuts-units-{nuts_year}.geojson")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule administrative_borders_nuts:
    message: "Normalise NUTS administrative borders."
    input:
        script = script_dir + "shapes/nuts.py",
        zipped = rules.download_raw_nuts_units_zip.output[0]
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


rule locations_template:
    message: "Generate locations configuration file for {wildcards.resolution} resolution from template."
    input:
        script = script_dir + "shapes/template_locations.py",
        template = model_template_dir + "locations.yaml",
        shapes = rules.units.output[0]
    output:
        yaml = "build/models/{resolution}/locations.yaml",
        csv = "build/models/{resolution}/locations.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/template_locations.py"


rule dwelling_density_shapes:
    message: "Create map of multi-family and single-family homes at the highest available spatial resolution"
    input:
        script = script_dir + "shapes/dwelling_map.py",
        dwellings = "data/automatic/eurostat-cens_11dwob_r3.tsv.gz",
        shapes = f"data/automatic/raw-nuts-units-{config['parameters']['eurostat-dataset-nuts-year']['cens_11dwob_r3']}.geojson",
    params:
        nuts_level = 3,
        crs = "EPSG:3035"
    conda: "../envs/geo.yaml"
    output: "build/data/shapes/dwellings-map.geojson"
    script: "../scripts/shapes/dwelling_map.py"


rule gross_added_value_shapes:
    message: "Create map of commercial sector gross added value at the highest available resolution"
    input:
        script = script_dir + "shapes/gva_map.py",
        eu_gva = "data/automatic/eurostat-nama_10r_3gva.tsv.gz",
        ch_gva = "build/data/ch-stats/gross_added_value.csv",
        shapes = f"data/automatic/raw-nuts-units-{config['parameters']['eurostat-dataset-nuts-year']['nama_10r_3gva']}.geojson",
    params:
        nuts_level = 3,
        crs = "EPSG:3035"
    conda: "../envs/geo.yaml"
    output: "build/data/shapes/gva-map.geojson"
    script: "../scripts/shapes/gva_map.py"


rule industry_employment_shapes:
    message: "Create map of {wildcards.industry_subsector} industry subsector employee numbers at the highest available resolution"
    input:
        script = script_dir + "shapes/employees_map.py",
        eu_business_statistics = "data/automatic/eurostat-sbs_r_nuts06_r2.tsv.gz",
        shapes = f"data/automatic/raw-nuts-units-{config['parameters']['eurostat-dataset-nuts-year']['sbs_r_nuts06_r2']}.geojson",
    params:
        nuts_level = 2,
        crs = "EPSG:3035"
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        industry_subsector = "FC_IND_MQ_E|FC_IND_FBT_E|FC_IND_NFM_E|FC_IND_TL_E|FC_IND_WP_E|FC_IND_PPP_E|FC_IND_CPC_E|FC_IND_NSP_E|FC_IND_NMM_E|FC_IND_IS_E|FC_IND_MAC_E|FC_IND_TE_E|FC_IND_CON_E"
    output: "build/data/shapes/{industry_subsector}-employees-map.geojson"
    script: "../scripts/shapes/employees_map.py"


rule industry_freight_shapes:
    message: "Create map of {wildcards.industry_subsector} industry subsector freight loading quantity at the highest available resolution"
    input:
        script = script_dir + "shapes/freight_map.py",
        eu_freight = "data/automatic/eurostat-road_go_na_rl3g.tsv.gz",
        shapes = f"data/automatic/raw-nuts-units-{config['parameters']['eurostat-dataset-nuts-year']['road_go_na_rl3g']}.geojson",
    params:
        nuts_level = 3,
        crs = "EPSG:3035"
    conda: "../envs/geo.yaml"
    wildcard_constraints:
        industry_subsector = "FC_IND_MQ_E|FC_IND_FBT_E|FC_IND_TL_E|FC_IND_WP_E|FC_IND_MAC_E|FC_IND_TE_E|FC_IND_NSP_E"
    output: "build/data/shapes/{industry_subsector}-freight-map.geojson"
    script: "../scripts/shapes/freight_map.py"


rule download_industry_sites_zip:
    message: "Download database pinpointing industry sites that are registered in the EU-ETS"
    params: url = config["data-sources"]["industry-sites"]
    output:
        protected("data/automatic/industry-sites.zip")
    conda: "../envs/shell.yaml"
    shell:
        "curl -sLo {output} '{params.url}'"


rule industry_sites:
    message: "Unzip industry sites."
    input: rules.download_industry_sites_zip.output[0]
    output: "build/data/shapes/Industrial_Database.csv"
    conda: "../envs/shell.yaml"
    shell: """
        unzip -j {input} "**/Industrial_Database.csv" -d build/data/shapes/
    """


rule industry_sites_points:
    message: "Geospatialise industry sites database"
    input:
        script = script_dir + "shapes/industry_emissions_map.py",
        industry_sites = rules.industry_sites.output[0]
    params:
        crs = "EPSG:3035"
    output: "build/data/shapes/{industry_subsector}-emissions-map.geojson"
    wildcard_constraints:
        industry_subsector = "FC_IND_IS_E|FC_IND_NMM_E|FC_IND_PPP_E|FC_IND_NFM_E|FC_IND_CPC_E|FC_IND_NSP_E"
    conda: "../envs/geo.yaml"
    script: "../scripts/shapes/industry_emissions_map.py"
