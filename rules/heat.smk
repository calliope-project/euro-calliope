rule download_gridded_weather_data:
    message: "Download gridded {wildcards.data_var} data"
    params: url = lambda wildcards: config["data-sources"]["gridded-weather-data"].format(data_var=wildcards.data_var)
    output: protected("data/automatic/gridded-weather/{data_var}.nc")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


rule download_when2heat_params:
    message: "Get parameters for heat demand profiles from the When2Heat project repository"
    output: directory("data/automatic/when2heat")
    params:
        url = lambda wildcards: config["data-sources"]["when2heat-params"].format(dataset=
            "{" + ",".join(["daily_demand.csv", "hourly_factors_COM.csv", "hourly_factors_MFH.csv", "hourly_factors_SFH.csv"]) + "}"
        )
    conda: "../envs/shell.yaml"
    shell: "mkdir -p {output} && curl -sSLo '{output}/#1' '{params.url}'"


rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        hh_end_use = "data/automatic/eurostat-hh-end-use.tsv.gz",
        ch_end_use = "data/automatic/ch-end-use.xlsx",
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "build/data/jrc-idees/tertiary/processed.csv",
        carrier_names = "config/energy-balances/energy-balance-carrier-names.csv"
    params:
        heat_tech_params = config["parameters"]["heat"],
        countries = config["scope"]["spatial"]["countries"],
        fill_missing_values = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
    conda: "../envs/default.yaml"
    output:
        total_demand = "build/data/heat/annual-heat-demand-twh.csv",
        electricity = "build/data/heat/annual-heat-electricity-demand-twh.csv",
    script: "../scripts/heat/annual_heat_demand.py"


rule rescale_annual_heat_demand_to_resolution:
    message: "Re-scale national heat demand at {wildcards.resolution} for household and commercial sectors"
    input:
        annual_demand = rules.annual_heat_demand.output["total_demand"],
        electricity = rules.annual_heat_demand.output["electricity"],
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv"
    conda: "../envs/default.yaml"
    output:
        total_demand = "build/data/heat/{resolution}/annual-heat-demand-twh.csv",
        electricity = "build/data/heat/{resolution}/annual-heat-electricity-demand-twh.csv",
    script: "../scripts/heat/rescale.py"


rule create_heat_demand_timeseries:
    message: "Create heat demand timeseries at {wildcards.resolution} for household and commercial sectors"
    input:
        annual_demand = rules.rescale_annual_heat_demand_to_resolution.output["total_demand"],
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        historic = False,
        power_scaling_factor = config["scaling-factors"]["power"],
    conda: "../envs/default.yaml"
    output:
        "build/models/{resolution}/timeseries/demand/electrified-heat-demand.csv",
    script: "../scripts/heat/create_timeseries.py"


use rule create_heat_demand_timeseries as create_heat_demand_timeseries_historic_electrification with:
    message: "Create timeseries for historic electrified heat demand"
    input:
        annual_demand = rules.rescale_annual_heat_demand_to_resolution.output["electricity"],
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        historic = True,
        power_scaling_factor = config["scaling-factors"]["power"],
    output:
        "build/models/{resolution}/timeseries/demand/heat-demand-historic-electrification.csv",


rule population_per_weather_gridbox:
    message: "Get {wildcards.resolution} population information per weather data gridbox"
    input:
        weather_grid = "data/automatic/gridded-weather/grid.nc",
        population = rules.raw_population_unzipped.output[0],
        locations = rules.units.output[0]
    params:
        lat_name = "lat",
        lon_name = "lon",
    conda: "../envs/geo.yaml"
    output: "build/data/{resolution}/population.nc"
    script: "../scripts/heat/population_per_gridbox.py"


rule unscaled_heat_profiles:
    message: "Generate gridded heat demand profile shapes for {wildcards.year} from weather and population data"

    input:
        population = rules.population_per_weather_gridbox.output[0],
        wind_speed = "data/automatic/gridded-weather/wind10m.nc",
        temperature = "data/automatic/gridded-weather/temperature.nc",
        when2heat = rules.download_when2heat_params.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
    conda: "../envs/default.yaml"
    output: "build/data/{resolution}/hourly_unscaled_heat_demand.nc"
    script: "../scripts/heat/unscaled_heat_profiles.py"
