"""Rules to process transport sector data."""


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        energy_balances="build/data/annual-energy-balances.csv",
        jrc_road_energy="build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance="build/data/jrc-idees/transport/processed-road-distance.csv",
        jrc_road_vehicles="build/data/jrc-idees/transport/processed-road-vehicles.csv",
    params:
        fill_missing_values=config["parameters"]["transport"]["fill-missing-values"],
        efficiency_quantile=config["parameters"]["transport"]["future-vehicle-efficiency-percentile"]
    conda: "../envs/default.yaml"
    output:
        distance="build/data/transport/annual-road-transport-distance-demand.csv",
        vehicles="build/data/transport/annual-road-transport-vehicles.csv",
        efficiency="build/data/transport/annual-road-transport-efficiency.csv",
        road_bau_electricity="build/data/transport/annual-road-transport-bau-electricity.csv",
    script: "../scripts/transport/annual_transport_demand.py"


rule create_road_transport_timeseries:
    message: "Create timeseries for road transport demand"
    input:
        data = "build/data/transport/annual-road-transport-distance-demand.csv",
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.type],
        type_name = lambda wildcards: config["parameters"]["transport"]["names"][wildcards.type],
        bau = False
    conda: "../envs/default.yaml"
    wildcard_constraints:
        type = "light-duty-vehicles|heavy-duty-vehicles|coaches-and-buses|passenger-cars|motorcycles"
    output:
        main = "build/data/transport/timeseries/timeseries-{type}.csv",
    script: "../scripts/transport/road_transport_timeseries.py"


use rule create_road_transport_timeseries as create_road_transport_timeseries_bau with:
    message: "Create timeseries for road transport demand"
    input:
        data = "build/data/transport/annual-road-transport-bau-electricity.csv"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        conversion_factor = lambda wildcards: config["parameters"]["transport"]["road-transport-conversion-factors"][wildcards.type],
        type_name = lambda wildcards: config["parameters"]["transport"]["names"][wildcards.type],
        bau = True
    output:
        "build/data/transport/timeseries/timeseries-{type}-bau.csv"


rule aggregate_timeseries: # TODO consider merge with other rules, as this is tiny atm
    message: "Aggregates timeseries for {wildcards.resolution} electrified road transport transport"
    input:
        time_series = (
            "build/data/transport/timeseries/timeseries-light-duty-vehicles.csv",
            "build/data/transport/timeseries/timeseries-heavy-duty-vehicles.csv",
            "build/data/transport/timeseries/timeseries-coaches-and-buses.csv",
            "build/data/transport/timeseries/timeseries-passenger-cars.csv",
            "build/data/transport/timeseries/timeseries-motorcycles.csv"),
    conda: "../envs/default.yaml"
    output:
        "build/models/{resolution}/timeseries/demand/electrified-road-transport.csv",
    script: "../scripts/transport/aggregate_timeseries.py"


use rule aggregate_timeseries as aggregate_bau_timeseries_bau with:
    message: "Aggregates timeseries for {wildcards.resolution} electrified road BAU transport"
    input:
        time_series=(
            "build/data/transport/timeseries/timeseries-light-duty-vehicles-bau.csv",
            "build/data/transport/timeseries/timeseries-coaches-and-buses-bau.csv",
            "build/data/transport/timeseries/timeseries-passenger-cars-bau.csv"),
    output:
        "build/models/{resolution}/timeseries/demand/electrified-bau-road-transport.csv"
