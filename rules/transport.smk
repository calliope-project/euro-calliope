"""Rules to process transport sector data."""

configfile: "./config/default.yaml"

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/transport/"


rule annual_road_transport_demand:
    message: "Calculate future road transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_road_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_road_energy = "build/data/jrc-idees/transport/processed-road-energy.csv",
        jrc_road_distance = "build/data/jrc-idees/transport/processed-road-distance.csv",
        jrc_road_vehicles = "build/data/jrc-idees/transport/processed-road-vehicles.csv",
    params:
        vehicle_efficiency_percentile = config["parameters"]["transport"]["future-vehicle-efficiency-percentile"]
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/data/transport/annual_road_transport_distance_demand.csv"),
        vehicles=temp("build/data/transport/annual_road_transport_vehicles.csv"),
        efficiency=temp("build/data/transport/annual_road_transport_efficiency.csv"),
        road_bau_electricity=temp("build/data/transport/annual_road_transport_bau_electricity.csv"),
    script: "../scripts/transport/annual_road_transport_demand.py"


rule annual_rail_transport_demand:
    message: "Calculate future rail transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_rail_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_rail_energy = "build/data/jrc-idees/transport/processed-rail-energy.csv",
        jrc_rail_distance = "build/data/jrc-idees/transport/processed-rail-distance.csv",
    conda: "../envs/default.yaml"
    output:
        rail_energy=temp("build/data/transport/annual_rail_transport_energy_demand.csv"),
        rail_bau_electricity=temp("build/data/transport/annual_rail_transport_bau_electricity.csv"),
    script: "../scripts/transport/annual_rail_transport_demand.py"


rule annual_air_transport_demand:
    message: "Calculate future air transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_air_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
    conda: "../envs/default.yaml"
    output:
        air_energy=temp("build/data/transport/annual_air_transport_energy_demand.csv"),
    script: "../scripts/transport/annual_air_transport_demand.py"


rule annual_marine_transport_demand:
    message: "Calculate future marine transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_marine_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
    conda: "../envs/default.yaml"
    output:
        marine_energy=temp("build/data/transport/annual_marine_transport_energy_demand.csv"),
    script: "../scripts/transport/annual_marine_transport_demand.py"
