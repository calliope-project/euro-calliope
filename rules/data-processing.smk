"""Rules to process input data."""

configfile: "./config/default.yaml"
localrules: eurostat_data_tsv, ch_data_xlsx
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule eurostat_data_tsv:
    message: "Get {wildcards.dataset} from Eurostat"
    params:
        url = lambda wildcards: config["data-sources"]["eurostat-base-url"].format(dataset=wildcards.dataset)
    output: protected("data/automatic/eurostat-{dataset}.tsv.gz")
    shell: "curl -sLo {output} {params.url}"


rule ch_data_xlsx:
    message: "Get {wildcards.dataset} from Swiss statistics"
    params:
        url = lambda wildcards: config["data-sources"]["swiss-stat"][wildcards.dataset]
    output: protected("data/automatic/ch-{dataset}.xlsx")
    shell: "curl -sLo {output} {params.url}"


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat and Switzerland-specific data"
    input:
        src = script_dir + "annual_energy_balance.py",
        eurostat_energy_balance = "data/automatic/eurostat-nrg_bal_c.tsv.gz",
        ch_energy_balance = "data/automatic/ch-energy-balance.xlsx",
        ch_industry_energy_balance = "data/automatic/ch-industry-energy-balance.xlsx",
        cat_names = config["data-sources"]["energy-balance-category-names"],
        carrier_names = config["data-sources"]["energy-balance-carrier-names"]
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    script: "../scripts/annual_energy_balance.py"


rule annual_road_transport_demand:
    message: "Calculate future road transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_road_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_road_energy = "data/transport/jrc_idees_processed_road_energy.csv",
        jrc_road_distance = "data/transport/jrc_idees_processed_road_distance.csv",
        jrc_road_vehicles = "data/transport/jrc_idees_processed_road_vehicles.csv",
    params:
        road_vehicle_efficiency = config["parameters"]["road-vehicle-efficiency"]
    conda: "../envs/default.yaml"
    output:
        distance=temp("build/annual_road_transport_distance_demand.csv"),
        vehicles=temp("build/annual_road_transport_vehicles.csv"),
        efficiency=temp("build/annual_road_transport_efficiency.csv"),
        road_bau_electricity=temp("build/annual_road_transport_bau_electricity.csv"),
    script: "../scripts/annual_road_transport_demand.py"


rule annual_rail_transport_demand:
    message: "Calculate future rail transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_rail_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
        jrc_rail_energy = "data/transport/jrc_idees_processed_rail_energy.csv",
        jrc_rail_distance = "data/transport/jrc_idees_processed_rail_distance.csv",
    conda: "../envs/default.yaml"
    output:
        rail_energy=temp("build/annual_rail_transport_energy_demand.csv"),
        rail_bau_electricity=temp("build/annual_rail_transport_bau_electricity.csv"),
    script: "../scripts/annual_rail_transport_demand.py"


rule annual_air_transport_demand:
    message: "Calculate future air transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_air_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
    conda: "../envs/default.yaml"
    output:
        air_energy=temp("build/annual_air_transport_energy_demand.csv"),  
    script: "../scripts/annual_air_transport_demand.py"


rule annual_marine_transport_demand:
    message: "Calculate future marine transport energy demand based on JRC IDEES"
    input:
        src = script_dir + "annual_marine_transport_demand.py",
        energy_balances = rules.annual_energy_balances.output[0],
    conda: "../envs/default.yaml"
    output:
        marine_energy=temp("build/annual_marine_transport_energy_demand.csv"),
    script: "../scripts/annual_marine_transport_demand.py"

