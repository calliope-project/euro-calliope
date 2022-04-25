"""Rules related to wind and solar."""

localrules: download_potentials, download_capacity_factors_wind_and_solar

ALL_WIND_AND_SOLAR_TECHNOLOGIES = [
    "wind-onshore", "wind-offshore", "open-field-pv",
    "rooftop-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat"
]


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


rule area_to_capacity_limits:
    message: "Use technology densities to convert wind & solar {wildcards.resolution} available area to capacity limits."
    input:
        script = script_dir + "wind-and-solar/capacity_limits.py",
        units = rules.units_without_shape.output[0],
        land_eligibility_km2 = rules.potentials.output.land_eligibility_km2,
    params:
        max_power_densities = config["parameters"]["maximum-installable-power-density"],
        roof_shares = config["parameters"]["roof-share"],
    output:
        rooftop = "build/data/{resolution}/supply/rooftop-solar.csv",
        offshore = "build/data/{resolution}/supply/wind-offshore.csv",
        onshore_and_open_field = "build/data/{resolution}/supply/open-field-solar-and-wind-onshore.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/wind-and-solar/capacity_limits.py"


rule capacity_factors_onshore_wind_and_solar:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        script = script_dir + "wind-and-solar/capacityfactors.py",
        locations = rules.units.output[0],
        timeseries = ancient("data/automatic/capacityfactors/{technology}-timeseries.nc"),
        coordinates = ancient("data/automatic/capacityfactors/wind-onshore-timeseries.nc")
    params:
        cf_threshold = config["capacity-factors"]["min"],
        gridcell_overlap_threshold=config["quality-control"]["capacity-factor-gridcell-overlap-threshold"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    wildcard_constraints:
        technology = "wind-onshore|rooftop-pv|open-field-pv|rooftop-pv-n|rooftop-pv-e-w|rooftop-pv-s-flat"
    output: "build/models/{resolution}/timeseries/supply/capacityfactors-{technology}.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/wind-and-solar/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for wind-offshore."
    input:
        script = script_dir + "wind-and-solar/capacityfactors_offshore.py",
        eez = rules.eez.output[0],
        shared_coast = rules.potentials.output.shared_coast,
        timeseries = ancient("data/automatic/capacityfactors/wind-offshore-timeseries.nc")
    params:
        cf_threshold = config["capacity-factors"]["min"],
        gridcell_overlap_threshold=config["quality-control"]["capacity-factor-gridcell-overlap-threshold"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    output: "build/models/{resolution}/timeseries/supply/capacityfactors-wind-offshore.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/wind-and-solar/capacityfactors_offshore.py"
