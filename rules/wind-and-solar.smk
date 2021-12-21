"""Rules related to wind and solar."""

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule download_potentials:
    message: "Download potential data."
    params: url = config["data-sources"]["potentials"]
    output: protected("data/automatic/raw-potentials.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule potentials:
    message: "Unzip potentials."
    input: rules.download_potentials.output[0]
    shadow: "minimal"
    output:
        land_eligibility_km2 = "build/data/{{resolution}}/{scenario}/areas.csv".format(
            scenario=config["parameters"]["wind-and-solar-potential-scenario"]
        ),
        shared_coast = "build/data/{resolution}/shared-coast.csv",
        demand = "build/data/{resolution}/demand.csv",
        population = "build/data/{resolution}/population.csv",
        land_cover = "build/data/{resolution}/land-cover.csv"
    conda: "../envs/shell.yaml"
    shell: "unzip -o {input} -d build/data"


rule download_capacity_factors_wind_and_solar:
    message: "Download data/automatic/capacityfactors/{wildcards.filename}."
    params: url = lambda wildcards: config["data-sources"]["capacity-factors"].format(filename=wildcards.filename)
    output: protected("data/automatic/capacityfactors/{filename}")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule capacity_factors_onshore_wind_and_solar:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        script = script_dir + "capacityfactors.py",
        locations = rules.units.output[0],
        timeseries = ancient("data/automatic/capacityfactors/{technology}-timeseries.nc"),
        coordinates = ancient("data/automatic/capacityfactors/wind-onshore-timeseries.nc")
    params:
        threshold = config["capacity-factors"]["min"],
        year = config["year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    wildcard_constraints:
        technology = "wind-onshore|rooftop-pv|open-field-pv|rooftop-pv-n|rooftop-pv-e-w|rooftop-pv-s-flat"
    output: "build/model/{resolution}/capacityfactors-{technology}.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for wind-offshore."
    input:
        script = script_dir + "capacityfactors_offshore.py",
        eez = rules.eez.output[0],
        shared_coast = rules.potentials.output.shared_coast,
        timeseries = ancient("data/automatic/capacityfactors/wind-offshore-timeseries.nc")
    params:
        threshold = config["capacity-factors"]["min"],
        year = config["year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    output: "build/model/{resolution}/capacityfactors-wind-offshore.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/capacityfactors_offshore.py"
