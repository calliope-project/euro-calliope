
rule annual_heat_demand:
    message: "Calculate national heat demand for household and commercial sectors"
    input:
        hh_end_use = "data/automatic/eurostat-hh-end-use.tsv.gz",
        ch_end_use = "data/automatic/ch-end-use.xlsx",
        energy_balance = rules.annual_energy_balances.output[0],
        commercial_demand = "build/data/jrc-idees/heat/commercial/jrc_idees_processed_energy.csv",
        carrier_names = "config/energy-balances/energy-balance-carrier-names.csv"
    params:
        heat_tech_params = config["parameters"]["heat-end-use"]
    conda: "../envs/default.yaml"
    output:
        demand=temp("build/data/heat/annual-heat-demand.csv"),
        electricity=temp("build/data/heat/annual-heat-electricity-consumption.csv"),
    script: "../scripts/heat/annual_heat_demand.py"

rule disagreggate_annual_heat_demand: # TODO have a subrule for the historic timeseries
    message: "Re-scale national heat demand at {wildcards.resolution} for household and commercial sectors"
    input:
        annual_demand = rules.annual_heat_demand.output["demand"],
        electricity = rules.annual_heat_demand.output["electricity"],
        populations = "build/data/regional/population.csv",
    params:
        locations = "config/locations/units_default_config_regional.csv", # FIXME we need all regions even in the minimal workflow, but this should be accessible outside pre-downloaded data
    conda: "../envs/default.yaml"
    output:
        demand = "build/data/heat/{resolution}/annual-heat-demand.csv",
        electricity = "build/data/heat/{resolution}/annual-heat-electricity-consumption.csv",
    script: "../scripts/heat/disaggregate_annual_heat_demand.py"


# rule create_heat_demand_timeseries: # TODO have realistic and separate space heat and water heat demand profiles
#     message: "Create heat demand timeseries for household and commercial sectors"
#     input:
#         annual_demand = rules.annual_heat_demand.output["demand"],
#         historic_electricity = rules.annual_heat_demand.output["electricity"],
#     params:
#         first_year = config["scope"]["temporal"]["first-year"],
#         final_year = config["scope"]["temporal"]["final-year"],
#         historic = False,
#         countries = config["scope"]["spatial"],
#         power_scaling_factor = config["scaling-factors"]["power"],
#     conda: "../envs/default.yaml"
#     wildcard_constraints
#     output:
#         main = "build/data/heat/timeseries/timeseries-heat-demand.csv",