"""Rules to process transport sector data."""


rule download_transport_timeseries:
    # TODO have correct timeseries data once RAMP has generated the new charging profile and it's been put on Zenodo
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


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        energy_balances = "build/data/annual-energy-balances.csv",
        jrc_road_energy = "build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.csv",
    params:
        fill_missing_values = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"],
        efficiency_quantile = config["parameters"]["transport"]["future-vehicle-efficiency-percentile"],
        uncontrolled_charging_share = config["parameters"]["transport"]["uncontrolled-ev-charging-share"],
    conda: "../envs/default.yaml"
    output:
        road_distance_controlled = "build/data/transport/annual-road-transport-distance-demand-controlled.csv",
        road_distance_uncontrolled = "build/data/transport/annual-road-transport-distance-demand-uncontrolled.csv",
        road_distance_historically_electrified = "build/data/transport/annual-road-transport-distance-demand-historic-electrification.csv",
    script: "../scripts/transport/annual_transport_demand.py"

rule create_controlled_road_transport_annual_demand_and_installed_capacities:
    message: "Create annual demand for controlled charging and corresponding charging potentials at {wildcards.resolution} resolution"
    input:
        annual_controlled_demand = "build/data/transport/annual-road-transport-distance-demand-controlled.csv",
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
        main = "build/data/{resolution}/demand/electrified-transport.csv",
    script: "../scripts/transport/road_transport_controlled_charging.py"

rule create_controlled_ev_charging_parameters:
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

rule create_uncontrolled_road_transport_timeseries:
    message: "Create timeseries for road transport demand  (uncontrolled charging)"
    input:
        annual_data = "build/data/transport/annual-road-transport-distance-demand-uncontrolled.csv",
        timeseries = "data/automatic/ramp-ev-consumption-profiles.csv.gz"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.vehicle_type],
        historic = False,
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
    conda: "../envs/default.yaml"
    wildcard_constraints:
        vehicle_type = "light-duty-vehicles|heavy-duty-vehicles|coaches-and-buses|passenger-cars|motorcycles"
    output:
        main = "build/data/transport/timeseries/timeseries-uncontrolled-{vehicle_type}.csv",
    script: "../scripts/transport/road_transport_timeseries.py"


use rule create_uncontrolled_road_transport_timeseries as create_uncontrolled_road_transport_timeseries_historic_electrification with:
    message: "Create timeseries for historic electrified road transport demand (uncontrolled charging)"
    input:
        annual_data = "build/data/transport/annual-road-transport-distance-demand-historic-electrification.csv",
        timeseries = "data/automatic/ramp-ev-consumption-profiles.csv.gz",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.vehicle_type],
        historic = True,
        countries = config["scope"]["spatial"]["countries"],
        country_neighbour_dict = config["data-pre-processing"]["fill-missing-values"]["ramp"],
    output:
        "build/data/transport/timeseries/timeseries-uncontrolled-{vehicle_type}-historic-electrification.csv"


rule aggregate_timeseries: # TODO consider merge with other rules, as this is tiny atm
    message: "Aggregates uncontrolled charging timeseries for {wildcards.resolution} electrified road transport transport"
    input:
        time_series = [
            f'build/data/transport/timeseries/timeseries-uncontrolled-{vehicle_type}.csv'
            for vehicle_type in config["parameters"]["transport"]["road-transport-conversion-factors"].keys()
        ],
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv"
    conda: "../envs/default.yaml"
    output:
        "build/models/{resolution}/timeseries/demand/uncontrolled-electrified-road-transport.csv",
    script: "../scripts/transport/aggregate_timeseries.py"


use rule aggregate_timeseries as aggregate_timeseries_historic_electrified with:
    message: "Aggregates uncontrolled charging timeseries for {wildcards.resolution} historically electrified road transport"
    input:
        time_series = (
            "build/data/transport/timeseries/timeseries-uncontrolled-light-duty-vehicles-historic-electrification.csv",
            "build/data/transport/timeseries/timeseries-uncontrolled-coaches-and-buses-historic-electrification.csv",
            "build/data/transport/timeseries/timeseries-uncontrolled-passenger-cars-historic-electrification.csv"),
        locations = "build/data/{resolution}/units.csv",
        populations = "build/data/{resolution}/population.csv"
    output:
        "build/models/{resolution}/timeseries/demand/uncontrolled-road-transport-historic-electrification.csv"
