"""Rules to process transport sector data."""


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        src=script_dir + "transport/annual_transport_demand.py",
        energy_balances="build/annual-energy-balances.csv",
        jrc_road_energy="build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance="build/data/jrc-idees/transport/processed-road-distance.csv",
        jrc_road_vehicles="build/data/jrc-idees/transport/processed-road-vehicles.csv",
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
        src=script_dir + "transport/road_transport_timeseries.py",
        road_distance_path="build/data/transport/annual-road-transport-distance-demand.csv",
        bau_electricity_path="build/data/transport/annual-road-transport-bau-electricity.csv"
    params:
        first_year=config["scope"]["temporal"]["first-year"],
        final_year=config["scope"]["temporal"]["final-year"],
        power_scaling_factor=config["scaling-factors"]["power"],
        ldv_conversion_factor=config["road_transport_conversion_factors"]["ldv_conversion_factor"],
        hdv_conversion_factor=config["road_transport_conversion_factors"]["hdv_conversion_factor"],
        coaches_and_buses_conversion_factor=config["road_transport_conversion_factors"][
            "coaches_and_buses_conversion_factor"],
        passenger_cars_conversion_factor=config["road_transport_conversion_factors"][
            "passenger_cars_conversion_factor"],
        powered_2_wheelers_conversion_factor=config["road_transport_conversion_factors"][
            "powered_2_wheelers_conversion_factor"]
    conda: "../envs/default.yaml"
    output:
        light_duty_vehicles_timeseries_out_path="build/data/transport/timeseries/timeseries-light-duty-vehicles.csv",
        heavy_duty_vehicles_timeseries_out_path="build/data/transport/timeseries/timeseries-heavy-duty-vehicles.csv",
        coaches_and_buses_timeseries_out_path="build/data/transport/timeseries/timeseries-coaches-and-buses.csv",
        passenger_cars_timeseries_out_path="build/data/transport/timeseries/timeseries-passenger-cars.csv",
        powered_2_wheelers_timeseries_out_path="build/data/transport/timeseries/timeseries-powered-2-wheelers.csv",
        light_duty_vehicles_bau_timeseries_out_path="build/data/transport/timeseries/timeseries-light-duty-vehicles-bau.csv",
        coaches_and_buses_bau_timeseries_out_path="build/data/transport/timeseries/timeseries-coaches-and-buses-bau.csv",
        passenger_cars_bau_timeseries_out_path="build/data/transport/timeseries/timeseries-passenger-cars-bau.csv",
    script: "../scripts/transport/road_transport_timeseries.py"


rule aggregate_timeseries:
    message: "Aggregates timeseries for {wildcards.resolution} electrified road transport and electrified road BAU transport"
    input:
        electrified_road_transport_timeseries=(
            "build/data/transport/timeseries/timeseries-light-duty-vehicles.csv",
            "build/data/transport/timeseries/timeseries-heavy-duty-vehicles.csv",
            "build/data/transport/timeseries/timeseries-coaches-and-buses.csv",
            "build/data/transport/timeseries/timeseries-passenger-cars.csv",
            "build/data/transport/timeseries/timeseries-powered-2-wheelers.csv"),
        electrified_road_bau_transport_timeseries=(
            "build/data/transport/timeseries/timeseries-light-duty-vehicles-bau.csv",
            "build/data/transport/timeseries/timeseries-coaches-and-buses-bau.csv",
            "build/data/transport/timeseries/timeseries-passenger-cars-bau.csv"),
    conda: "../envs/default.yaml"
    output:
        road_transport_timeseries="build/models/{resolution}/timeseries/demand/electrified-road-transport.csv",
        road_transport_bau_timeseries="build/models/{resolution}/timeseries/demand/electrified-bau-road-transport.csv"
    script: "../scripts/transport/aggregate_timeseries.py"
