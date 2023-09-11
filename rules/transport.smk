"""Rules to process transport sector data."""


rule annual_transport_demand:
    message: "Calculate future transport energy demand based on JRC IDEES"
    input:
        src = "src/construct/annual_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_road_energy = "data/transport/jrc_idees_processed_road_energy.csv",
        jrc_road_distance = "data/transport/jrc_idees_processed_road_distance.csv",
        jrc_road_vehicles = "data/transport/jrc_idees_processed_road_vehicles.csv",
        jrc_rail_energy = "data/transport/jrc_idees_processed_rail_energy.csv",
        jrc_rail_distance = "data/transport/jrc_idees_processed_rail_distance.csv",
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/annual_road_transport_distance_demand.csv"),
        vehicles=temp("build/annual_road_transport_vehicles.csv"),
        efficiency=temp("build/annual_road_transport_efficiency.csv"),
        rail_energy=temp("build/annual_rail_transport_energy_demand.csv"),
        air_energy=temp("build/annual_air_transport_energy_demand.csv"),
        marine_energy=temp("build/annual_marine_transport_energy_demand.csv"),
        road_bau_electricity=temp("build/annual_road_transport_bau_electricity.csv"),
        rail_bau_electricity=temp("build/annual_rail_transport_bau_electricity.csv"),
    script: "../src/construct/annual_transport_demand.py"
