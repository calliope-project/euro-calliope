"""Rules to process transport sector data."""


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "transport/annual_transport_demand.py",
        energy_balances = "build/annual_energy_balances.csv",
        jrc_road_energy = "build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.csv",
        jrc_road_vehicles = "build/data/jrc-idees/transport/processed-road-vehicles.csv",
    conda: "../envs/default.yaml"
    output:
        distance="build/data/new/annual_road_transport_distance_demand.csv",
        vehicles="build/data/new/annual_road_transport_vehicles.csv",
        efficiency="build/data/new/annual_road_transport_efficiency.csv",
        road_bau_electricity="build/data/new/annual_road_transport_bau_electricity.csv",
    script: "../scripts/transport/annual_transport_demand.py"


rule create_road_transport_demand_timeseries:
    message: "Create timeseries for road transport demand"
    input:
        src = script_dir + "transport/road_transport_timeseries.py",
        road_distance_path = "build/data/new/annual_road_transport_distance_demand.csv",
        bau_electricity_path = "build/data/new/annual_road_transport_bau_electricity.csv"
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        power_scaling_factor = config["scaling-factors"]["power"],
        ldv_conversion_factor = config["road_transport_conversion_factors"]["ldv_conversion_factor"],
        hdv_conversion_factor = config["road_transport_conversion_factors"]["hdv_conversion_factor"],
        coaches_and_buses_conversion_factor = config["road_transport_conversion_factors"]["coaches_and_buses_conversion_factor"],
        passenger_cars_conversion_factor = config["road_transport_conversion_factors"]["passenger_cars_conversion_factor"],
        powered_2_wheelers_conversion_factor = config["road_transport_conversion_factors"]["powered_2_wheelers_conversion_factor"]
    conda: "../envs/default.yaml"
    output:
        light_duty_vehicles_timeseries_out_path="build/data/new/timeseries_light_duty_vehicles.csv",
        heavy_duty_vehicles_timeseries_out_path="build/data/new/timeseries_heavy_duty_vehicles.csv",
        coaches_and_buses_timeseries_out_path="build/data/new/timeseries_coaches_and_buses.csv",
        passenger_cars_timeseries_out_path="build/data/new/timeseries_passenger_cars.csv",
        powered_2_wheelers_timeseries_out_path="build/data/new/timeseries_powered_2_wheelers.csv",
        light_duty_vehicles_bau_timeseries_out_path="build/data/new/timeseries_light_duty_vehicles_bau.csv",
        coaches_and_buses_bau_timeseries_out_path="build/data/new/timeseries_coaches_and_buses_bau.csv",
        passenger_cars_bau_timeseries_out_path="build/data/new/timeseries_passenger_cars_bau.csv",


    script: "../scripts/transport/road_transport_timeseries.py"
