"""Rules to generate hydro electricity capacities and time series."""

configfile: "./config/default.yaml"
localrules: download_runoff_data, download_stations_database, stations_database
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule download_runoff_data:
    message: "Create an atlite cutout of Europe consisting of ERA5 runoff data."
    input:
        script = script_dir + "hydro/runoff.py"
    params:
        year = config["year"],
        x_min = config["scope"]["bounds"]["x_min"],
        x_max = config["scope"]["bounds"]["x_max"],
        y_min = config["scope"]["bounds"]["y_min"],
        y_max = config["scope"]["bounds"]["y_max"]
    output:
        protected("data/automatic/europe-cutout-{}.nc".format(config["year"]))
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/runoff.py"


rule download_stations_database:
    message: "Download database of hydro electricity stations."
    params: url = config["data-sources"]["hydro-stations"]
    output:
        protected("data/automatic/raw-hydro-stations.zip")
    shell:
        "curl -sLo {output} '{params.url}'"


rule stations_database:
    message: "Unzip stations database."
    input: rules.download_stations_database.output
    output: "build/data/jrc-hydro-power-plant-database.csv"
    shadow: "full"
    shell:
        """
        unzip {input} -d ./build/
        mv build/energy-modelling-toolkit-hydro-power-database-f616a8d/data/jrc-hydro-power-plant-database.csv {output}
        """


rule fix_basins:
    message: "Fix invalid basins."
    input:
        script = script_dir + "hydro/fix_basins.py",
        basins = config["data-sources"]["hydro-basins"]
    output: "build/data/hybas_eu_lev07_v1c.gpkg"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/fix_basins.py"


rule preprocess_hydro_stations:
    # Some locations of stations are imprecise and in the sea. Slightly move them.
    # Some other stations seem incorrect. Remove.
    # Add missing pumped hydro stations in Romania.
    message: "Preprocess hydro stations."
    input:
        script = script_dir + "hydro/preprocess_hydro_stations.py",
        stations = rules.stations_database.output[0],
        basins = rules.fix_basins.output[0]
    params: buffer_size = 1 / 60 # move stations up to 1 arcminute < 1 km
    output: "build/data/jrc-hydro-power-plant-database-preprocessed.csv"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/preprocess_hydro_stations.py"


rule inflow_m3:
    message: "Determine water inflow time series for all hydro electricity."
    input:
        script = script_dir + "hydro/inflow_m3.py",
        stations = rules.preprocess_hydro_stations.output[0],
        basins = rules.fix_basins.output[0],
        runoff = rules.download_runoff_data.output[0]
    params: year = config["year"]
    output: "build/data/hydro-electricity-with-water-inflow.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_m3.py"


rule inflow_mwh:
    message: "Determine energy inflow time series for all hydro electricity."
    input:
        script = script_dir + "hydro/inflow_mwh.py",
        stations = rules.inflow_m3.output[0],
        generation = config["data-sources"]["irena-generation"]
    params:
        year = config["year"],
        max_capacity_factor = config["capacity-factors"]["max"]
    output: "build/data/hydro-electricity-with-energy-inflow.nc"
    conda: "../envs/hydro.yaml"
    script: "../scripts/hydro/inflow_mwh.py"
