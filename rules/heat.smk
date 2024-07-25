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

rule download_heat_pump_characteristics:
    message: "Download manufacturer heat pump data"
    params: url = config["data-sources"]["heat-pump-characteristics"]
    output: protected("data/automatic/heat-pump-characteristics.nc")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        hh_end_use = "data/automatic/eurostat-hh-end-use.tsv.gz",
        ch_end_use = "data/automatic/ch-end-use.xlsx",
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "build/data/jrc-idees/tertiary/processed.csv",
        carrier_names = "config/energy-balances/energy-balance-carrier-names.csv"
    params:
        heat_tech_params = config["parameters"]["heat"]["tech-efficiencies"],
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
        electricity = "build/data/heat/{resolution}/annual-historic-electrified-heat-demand-twh.csv",
    script: "../scripts/heat/rescale.py"


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
    message: "Generate gridded heat demand profile shapes from weather and population data"
    input:
        wind_speed = "data/automatic/gridded-weather/wind10m.nc",
        temperature = "data/automatic/gridded-weather/temperature.nc",
        when2heat = rules.download_when2heat_params.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
    conda: "../envs/default.yaml"
    output: "build/data/heat/hourly_unscaled_heat_demand.nc"
    script: "../scripts/heat/unscaled_heat_profiles.py"


rule heat_pump_cop:
    message: "Generate gridded heat pump coefficient of performance (COP)"
    input:
        temperature_air = "data/automatic/gridded-weather/temperature.nc",
        temperature_ground = "data/automatic/gridded-weather/tsoil5.nc",
        heat_pump_characteristics = rules.download_heat_pump_characteristics.output[0]
    params:
        sink_temperature = config["parameters"]["heat-pump"]["sink-temperature"],
        space_heat_sink_shares = config["parameters"]["heat-pump"]["space-heat-sink-shares"],
        correction_factor = config["parameters"]["heat-pump"]["correction-factor"],
        heat_pump_shares = config["parameters"]["heat-pump"]["heat-pump-shares"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
    conda: "../envs/default.yaml"
    output: "build/data/heat/heat-pump-cop.nc"
    script: "../scripts/heat/heat_pump_cop.py"

rule group_gridded_timeseries:
    message: "Generate {wildcards.resolution} {wildcards.input_dataset} timeseries data from gridded data "
    input:
        gridded_timeseries_data = "build/data/heat/{input_dataset}.nc",
        grid_weights = rules.population_per_weather_gridbox.output[0],
    conda: "../envs/default.yaml"
    threads: 4
    output: temp("build/data/heat/{resolution}/{input_dataset}.nc")
    script: "../scripts/heat/group_gridded_timeseries.py"


rule heat_pump_final_timeseries:
    message: "Combine hot water and space heating characteristics to generate a weighted average {wildcards.resolution} {wildcards.input_dataset} `heat` carrier timeseries."
    input:
        timeseries_data = "build/data/heat/{resolution}/{input_dataset}.nc",
        annual_demand = rules.rescale_annual_heat_demand_to_resolution.output.total_demand
    conda: "../envs/default.yaml"
    wildcard_constraints:
        input_dataset = "heat-pump-cop"
    output: "build/models/{resolution}/timeseries/conversion/{input_dataset}.csv"
    script: "../scripts/heat/heat_pump_final_timeseries.py"


rule heat_demand_final_timeseries:
    message: "Generate {wildcards.resolution} {wildcards.input_dataset} {wildcards.tech_group} timeseries data from gridded data "
    input:
        timeseries_data = "build/data/heat/{resolution}/hourly_unscaled_heat_demand.nc",
        annual_demand = "build/data/heat/{resolution}/annual-{input_dataset}-demand-twh.csv",
    conda: "../envs/default.yaml"
    params:
        sfh_mfh_shares = config["parameters"]["heat"]["sfh-mfh-shares"],
        scaling_factor = config["scaling-factors"]["power"]
    wildcard_constraints:
        input_dataset = "heat|historic-electrified-heat",
        tech_group = "demand|supply"
    output: "build/models/{resolution}/timeseries/{tech_group}/{input_dataset}.csv"
    script: "../scripts/heat/heat_demand_final_timeseries.py"

use rule heat_demand_final_timeseries as electrified_heat_demand_final_timeseries with:
    input:
        timeseries_data = "build/data/heat/{resolution}/hourly_unscaled_heat_demand.nc",
        annual_demand = "build/data/heat/{resolution}/annual-heat-demand-twh.csv",
        cop = "build/data/heat/{resolution}/heat-pump-cop.nc"
    params:
        sfh_mfh_shares = config["parameters"]["heat"]["sfh-mfh-shares"],
        scaling_factor = config["scaling-factors"]["power"],
        electrification_shares = config["parameters"]["heat"]["electrification-shares"]
    wildcard_constraints:
        input_dataset = "electrified-heat",
        tech_group = "demand"
