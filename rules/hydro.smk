"""Rules to generate hydro electricity capacities and time series."""

localrules: download_basins_database, download_stations_database
localrules: download_hydro_generation_data, download_pumped_hydro_data, basins_database, stations_database


rule download_hydro_generation_data:
    message: "Download database of historical hydro power generation."
    params: url = config["data-sources"]["hydro-generation"]
    output:
        protected("data/automatic/raw-hydro-generation.csv")
    conda: "../envs/shell.yaml"
    shell:
        "curl -sLo {output} '{params.url}'"


rule download_pumped_hydro_data:
    message: "Download database of pumped hydro storage capacity data."
    params: url = config["data-sources"]["national-phs-storage-capacities"]
    output:
        protected("data/automatic/raw-pumped-hydro-storage-capacities-gwh.csv")
    conda: "../envs/shell.yaml"
    shell:
        "curl -sLo {output} '{params.url}'"


rule download_runoff_data:
    message: "Create an atlite cutout of Europe consisting of ERA5 runoff data between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        script = script_dir + "hydro/runoff.py"
    params:
        x_min = config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = config["scope"]["spatial"]["bounds"]["y_max"]
    output:
        protected("data/automatic/europe-cutout-{first_year}-{final_year}.nc")
    conda: "../envs/hydro.yaml"
    envmodules: "eth_proxy"
    script: "../scripts/hydro/runoff.py"


rule download_basins_database:
    message: "Download database of hydro basins."
    params: url = config["data-sources"]["hydro-basins"]
    output:
        protected("data/automatic/raw-hydro-basins.zip")
    conda: "../envs/shell.yaml"
    shell:
        "curl -sLo {output} '{params.url}'"


rule download_stations_database:
    message: "Download database of hydro electricity stations."
    params: url = config["data-sources"]["hydro-stations"]
    output:
        protected("data/automatic/raw-hydro-stations.zip")
    conda: "../envs/shell.yaml"
    shell:
        "curl -sLo {output} '{params.url}'"


rule basins_database:
    message: "Unzip basins database."
    input: rules.download_basins_database.output
    output: "build/data/basins/hybas_eu_lev07_v1c.shp"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} -d ./build/data/basins/"


rule stations_database:
    message: "Unzip stations database."
    input: rules.download_stations_database.output
    output: "build/data/jrc-hydro-power-plant-database.csv"
    shadow: "full"
    conda: "../envs/shell.yaml"
    shell:
        """
        unzip -j {input} "**/jrc-hydro-power-plant-database.csv" -d build/data/
        """


rule preprocess_basins:
    message: "Preprocess basins."
    input:
        script = script_dir + "hydro/preprocess_basins.py",
        basins = rules.basins_database.output[0]
    params:
        x_min = config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = config["scope"]["spatial"]["bounds"]["y_max"]
    output: "build/data/hybas_eu_lev07_v1c.gpkg"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/preprocess_basins.py"


rule preprocess_hydro_stations:
    message: "Preprocess hydro stations."
    input:
        script = script_dir + "hydro/preprocess_hydro_stations.py",
        stations = rules.stations_database.output[0],
        basins = rules.preprocess_basins.output[0],
        phs_storage_capacities = rules.download_pumped_hydro_data.output[0]
    params:
        buffer_size_m = config["quality-control"]["hydro"]["station-nearest-basin-max-km"] * 1000,
        countries = config["scope"]["spatial"]["countries"],
        scale_phs = config["quality-control"]["hydro"]["scale-phs-according-to-geth-et-al"]
    output: "build/data/jrc-hydro-power-plant-database-preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/preprocess_hydro_stations.py"


rule inflow_m3:
    message: "Determine water inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        script = script_dir + "hydro/inflow_m3.py",
        stations = rules.preprocess_hydro_stations.output[0],
        basins = rules.preprocess_basins.output[0],
        runoff = rules.download_runoff_data.output[0]
    output: "build/data/hydro-electricity-with-water-inflow-{first_year}-{final_year}.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for all hydro electricity between the years {wildcards.first_year} and {wildcards.final_year}."
    input:
        script = script_dir + "hydro/inflow_mwh.py",
        stations = rules.inflow_m3.output[0],
        generation = rules.download_hydro_generation_data.output[0]
    params:
        max_capacity_factor = config["capacity-factors"]["max"]
    output: "build/data/hydro-electricity-with-energy-inflow-{first_year}-{final_year}.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_mwh.py"


rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        script = script_dir + "hydro/hydro_capacities.py",
        locations = rules.units.output[0],
        plants = rules.preprocess_hydro_stations.output[0]
    output: "build/data/{resolution}/hydro-capacities-mw.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/hydro/hydro_capacities.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        script = script_dir + "hydro/capacityfactors_hydro.py",
        capacities = rules.hydro_capacities.output[0],
        stations = "build/data/hydro-electricity-with-energy-inflow-{first_year}-{final_year}.nc".format(
            first_year = config["scope"]["temporal"]["first-year"],
            final_year = config["scope"]["temporal"]["final-year"]
        ),
        locations = rules.units.output[0]
    params:
        threshold = config["capacity-factors"]["min"]
    output:
        ror = "build/models/{resolution}/timeseries/supply/capacityfactors-hydro-ror.csv",
        reservoir = "build/models/{resolution}/timeseries/supply/capacityfactors-hydro-reservoir-inflow.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/hydro/capacityfactors_hydro.py"


rule hydro_storage_techs_at_locations_template:
    message: "Allocate hydro {wildcards.hydro_tech_type} techs to {wildcards.resolution} locations from template."
    input:
        script = script_dir + "template_techs.py",
        template = techs_template_dir + "{hydro_tech_type}/hydro.yaml",
        locations = rules.hydro_capacities.output[0],
        timeseries_data = rules.capacity_factors_hydro.output,
    params:
        capacity_factors = config["capacity-factors"]["average"],
        scaling_factors = config["scaling-factors"],
    wildcard_constraints:
        hydro_tech_type = "storage|supply"
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/{hydro_tech_type}/hydro.yaml"
    script: "../scripts/template_techs.py"
