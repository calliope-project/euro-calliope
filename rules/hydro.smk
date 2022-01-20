"""Rules to generate hydro electricity capacities and time series."""

localrules: download_basins_database, download_stations_database
localrules: download_hydro_generation_data, download_pumped_hydro_data, basins_database, stations_database
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


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
    message: "Create an atlite cutout of Europe consisting of ERA5 runoff data."
    input:
        script = script_dir + "hydro/runoff.py"
    params:
        # Need an extra initial year of data, since runoff inflow is shifted in time by atlite in `inflow_m3`
        first_year = config["scope"]["temporal"]["first-year"] - 1,
        final_year = config["scope"]["temporal"]["final-year"],
        x_min = config["scope"]["spatial"]["bounds"]["x_min"],
        x_max = config["scope"]["spatial"]["bounds"]["x_max"],
        y_min = config["scope"]["spatial"]["bounds"]["y_min"],
        y_max = config["scope"]["spatial"]["bounds"]["y_max"]
    output:
        protected("data/automatic/europe-cutout-{params.first_year}-{params.final_year}.nc")
    conda: "../envs/hydro.yaml"
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
    message: "Determine water inflow time series for all hydro electricity."
    input:
        script = script_dir + "hydro/inflow_m3.py",
        stations = rules.preprocess_hydro_stations.output[0],
        basins = rules.preprocess_basins.output[0],
        runoff = rules.download_runoff_data.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
    output: "build/data/hydro-electricity-with-water-inflow.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for all hydro electricity."
    input:
        script = script_dir + "hydro/inflow_mwh.py",
        stations = rules.inflow_m3.output[0],
        generation = rules.download_hydro_generation_data.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        max_capacity_factor = config["capacity-factors"]["max"]
    output: "build/data/hydro-electricity-with-energy-inflow.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_mwh.py"
