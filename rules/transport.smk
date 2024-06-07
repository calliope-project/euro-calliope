"""Rules to process transport sector data."""


rule download_transport_timeseries:
    message: "Get EV data from RAMP"
    params:
        url = config["data-sources"]["controlled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output:
        protected("data/automatic/ramp-ev-{dataset}.csv.gz")
    wildcard_constraints:
        dataset = "consumption-profiles|plugin-profiles"
    localrule: True
    shell: "curl -sSLo {output} {params.url}"

rule download_uncontrolled_timeseries:
    # TODO: move into rule download_transport_timeseries once PR 356 is merged
    message: "Get EV uncontrolled charging data from RAMP"
    params:
        url = config["data-sources"]["uncontrolled-ev-profiles"]
    conda: "../envs/shell.yaml"
    output: protected("data/automatic/ramp-ev-uncontrolled-charging-profiles.csv.gz")
    localrule: True
    shell: "curl -sSLo {output} {params.url}"

rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        energy_balances = "build/data/annual-energy-balances.csv",
        jrc_road_energy = "build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.csv",
    params:
        fill_missing_values = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"],
    output:
        road_distance = "build/data/transport/annual-road-transport-distance-demand.csv",
        road_distance_historically_electrified = "build/data/transport/annual-road-transport-distance-demand-historic-electrification.csv",
    conda: "../envs/default.yaml"
    script: "../scripts/transport/annual_transport_demand.py"

rule create_road_transport_vehicle_parameters:
    message: "Create vehicle parameters at {wildcards.resolution} resolution"
    input:
        ev_vehicle_number = "build/data/jrc-idees/transport/processed-road-vehicles.csv",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.csv",
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        transport_scaling_factor = config["scaling-factors"]["transport"],
        battery_sizes = config["parameters"]["transport"]["ev-battery-sizes"],
        conversion_factors = config["parameters"]["transport"]["road-transport-conversion-factors"],
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
    conda: "../envs/default.yaml"
    output:
        main = "build/data/{resolution}/supply/transport.csv",
    script: "../scripts/transport/road_transport_vehicle_parameters.py"

rule create_transport_demand_data_for_yaml:
    message: "Create transport demand data to be exported into the YAML file"
    input:
        annual_road_distance = "build/data/transport/annual-road-transport-distance-demand.csv",
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        countries = config["scope"]["spatial"]["countries"],
        vehicle_aggregation = config["parameters"]["transport"]["vehicle-type-aggregation"],
        transport_factor = config["scaling-factors"]["transport"],
    conda: "../envs/default.yaml"
    output:
        main = "build/data/{resolution}/demand/transport.csv",
    script: "../scripts/transport/road_transport_demand_data_for_yaml.py"

rule create_controlled_ev_charging_parameters:
    # Necessary timeseries for custom constraints of controlled charging
    message: "Create timeseries parameters {wildcards.dataset_name} for controlled EV charging at {wildcards.resolution} resolution"
    input:
        ev_profiles = lambda wildcards: "data/automatic/ramp-ev-consumption-profiles.csv.gz" if "demand" in wildcards.dataset_name else f"data/automatic/ramp-ev-{wildcards.dataset_name}.csv.gz",
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv",
    params:
        demand_range = config["parameters"]["transport"]["monthly-demand-bound-fraction"],
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
        countries = config["scope"]["spatial"]["countries"],
    wildcard_constraints:
        dataset_name = "demand-shape-equals|demand-shape-max|demand-shape-min|plugin-profiles"
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/timeseries/demand/{dataset_name}-ev.csv"
    script: "../scripts/transport/road_transport_controlled_constraints.py"

rule create_uncontrolled_road_transport_timeseries_historic_electrification:
    message: "Create timeseries for historical electrified road transport demand  (assumed uncontrolled charging)"
    input:
        annual_data = "build/data/transport/annual-road-transport-distance-demand-historic-electrification.csv",
        timeseries = "data/automatic/ramp-ev-uncontrolled-charging-profiles.csv.gz"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.vehicle_type],
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
    conda: "../envs/default.yaml"
    wildcard_constraints:
        vehicle_type = "light-duty-vehicles|heavy-duty-vehicles|coaches-and-buses|passenger-cars|motorcycles"
    output:
        main = "build/data/transport/timeseries/timeseries-uncontrolled-{vehicle_type}-historic-electrification.csv",
    script: "../scripts/transport/road_transport_historic_electrification.py"


rule aggregate_timeseries_historic_electrified: # TODO consider merge with other rules, as this is tiny atm
    message: "Aggregates uncontrolled charging timeseries for {wildcards.resolution} electrified road transport transport"
    input:
        time_series = (  # Historical data does not have electrified motorcycles and heavy-duty-vehicles
            "build/data/transport/timeseries/timeseries-uncontrolled-light-duty-vehicles-historic-electrification.csv",
            "build/data/transport/timeseries/timeseries-uncontrolled-coaches-and-buses-historic-electrification.csv",
            "build/data/transport/timeseries/timeseries-uncontrolled-passenger-cars-historic-electrification.csv",
        ),
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv"
    conda: "../envs/default.yaml"
    output:
        "build/models/{resolution}/timeseries/demand/uncontrolled-road-transport-historic-electrification.csv",
    script: "../scripts/transport/aggregate_timeseries.py"