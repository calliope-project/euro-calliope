"""Rules to process transport sector data."""


rule annual_road_transport_demand:
    message: "Calculate future road transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "transport/annual_road_rail_transport_demand.py",
        energy_balances = "build/data/eurostat/annual-road-transport-energy-balances.nc",
        jrc_road_energy = "build/data/jrc-idees/transport/processed-road-energy.nc",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.nc",
        jrc_road_vehicles = "build/data/jrc-idees/transport/processed-road-vehicles.nc",
    params:
        vehicle_efficiency_percentile = config["parameters"]["transport"]["future-vehicle-efficiency-percentile"]
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/data/transport/annual_road_transport_distance_demand.csv"),
        vehicles=temp("build/data/transport/annual_road_transport_vehicles.csv"),
        efficiency=temp("build/data/transport/annual_road_transport_efficiency.csv"),
        road_bau_electricity=temp("build/data/transport/annual_road_transport_bau_electricity.csv"),
    script: "../scripts/transport/annual_road_rail_transport_demand.py"
