"""Rules to process input data."""

configfile: "./config/default.yaml"
localrules: eurostat_data_tsv, ch_data_xlsx
root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule eurostat_data_tsv:
    message: "Get various datasets from Eurostat"
    output:
        energy_balance = protected("data/automatic/annual_energy_balances.tsv.gz"),
        hh_end_use = protected("data/automatic/hh_end_use.tsv.gz"),
        freight = protected("data/automatic/freight.tsv.gz"),
        employees = protected("data/automatic/employees.tsv.gz"),
        gva = protected("data/automatic/gva.tsv.gz"),
        dwellings = protected("data/automatic/dwellings.tsv.gz"),
    shell:
        """
        curl -sLo {output.energy_balance} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][energy-balance]}.tsv.gz'
        curl -sLo {output.hh_end_use} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][hh-end-use]}.tsv.gz'
        curl -sLo {output.freight} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][freight]}.tsv.gz'
        curl -sLo {output.employees} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][employees]}.tsv.gz'
        curl -sLo {output.gva} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][gva]}.tsv.gz'
        curl -sLo {output.dwellings} '{config[data-sources][eurostat][base-url]}{config[data-sources][eurostat][dwellings]}.tsv.gz'
        """


rule ch_data_xlsx:
    message: "Get Swiss annual energy balances and household end uses"
    output:
        energy_balance = protected("data/automatic/ch_annual_energy_balances.xlsx"),
        industry_energy_balance = protected("data/automatic/ch_annual_industry_energy_balances.xlsx"),
        end_use = protected("data/automatic/ch_hh_end_use.xlsx"),
        gva = protected("data/automatic/ch_gva.xlsx")
    shell:
        """
        curl -sLo {output.energy_balance} '{config[data-sources][swiss-stat][energy-balance]}'
        curl -sLo {output.industry_energy_balance} '{config[data-sources][swiss-stat][industry-energy-balance]}'
        curl -sLo {output.end_use} '{config[data-sources][swiss-stat][end-use]}'
        curl -sLo {output.gva} '{config[data-sources][swiss-stat][gva]}'
        """


rule annual_energy_balances:
    message: "Process annual energy balances from Eurostat and Switzerland-specific data"
    input:
        src = script_dir + "annual_energy_balance.py",
        energy_balance = rules.eurostat_data_tsv.output.energy_balance,
        ch_energy_balance = rules.ch_data_xlsx.output.energy_balance,
        ch_industry_energy_balance = rules.ch_data_xlsx.output.industry_energy_balance,
        cat_names = "data/energy_balance_category_names.csv",
        carrier_names = "data/energy_balance_carrier_names.csv"
    output: "build/annual_energy_balances.csv"
    params:
        countries = config["scope"]["countries"]
    conda: "../envs/default.yaml"
    shadow: "minimal"
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

