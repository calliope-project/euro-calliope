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
