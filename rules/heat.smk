rule commercial_heat_supply_by_end_use:
    message: "Assign Eurostat annual {wildcards.building_sector} fuel supply to building heat end uses using distributions from JRC-IDEES"
    input:
        script = script_dir + "heat/annual_commercial_industry_heat_energy_supply.py",
        ch_building_heat_energy_supply="build/data/ch-stats/building-heat-energy-balance.nc",
        annual_energy_balances="build/data/eurostat/annual-{building_sector}-energy-balances.nc",
        jrc_idees_sector_energy_supply="build/data/jrc-idees/tertiary/processed-energy.nc"
    params:
        countries = config["scope"]["spatial"]["countries"]
    wildcard_constraints:
        building_sector = "commercial"
    conda: "../envs/default.yaml"
    output: "build/data/heat/annual-{building_sector}-heat-energy-supply.nc"
    script: "../scripts/heat/annual_commercial_industry_heat_energy_supply.py"


use rule commercial_heat_supply_by_end_use as industry_heat_supply_by_end_use with:
    input:
        script = script_dir + "heat/annual_commercial_industry_heat_energy_supply.py",
        annual_energy_balances = "build/data/eurostat/annual-{building_sector}-energy-balances.nc",
        jrc_idees_sector_energy_supply = "build/data/jrc-idees/{building_sector}/processed-energy.nc"
    wildcard_constraints:
        building_sector = "industry"


rule household_heat_supply_by_end_use:
    message: "Assign Eurostat annual commercial fuel supply to building heat end uses using distributions from JRC-IDEES"
    input:
        script = script_dir + "heat/annual_household_heat_energy_supply.py",
        ch_building_heat_energy_supply="build/data/ch-stats/building-heat-energy-balance.nc",
        annual_energy_balances="build/data/eurostat/annual-household-energy-balances.nc",
        eurostat_household_building_heat="build/data/eurostat/household-building-heat-end-use-energy-balances.nc"
    params:
        countries = config["scope"]["spatial"]["countries"],
        carrier_names = config["mapping-keys"]["eurostat"]["carrier-names"]
    conda: "../envs/default.yaml"
    output: "build/data/heat/annual-household-heat-energy-supply.nc"
    script: "../scripts/heat/annual_household_heat_energy_supply.py"


rule building_heat_demand:
    message: "Derive {wildcards.demand_type} from annual energy consumption to meet heat demand in {wildcards.building_sector} buildings."
    input:
        script = script_dir + "heat/annual_heat_demand.py",
        annual_energy_balances = "build/data/eurostat/annual-{building_sector}-energy-balances.nc",
        end_use_energy_supply = "build/data/heat/annual-{building_sector}-heat-energy-supply.nc",
    conda: "../envs/default.yaml"
    params:
        heat_tech_efficiency_overrides = config["parameters"]["heat"]["supply-to-demand-efficiencies"]
    wildcard_constraints:
        demand_type = "heat-demand|current-electricity-consumption"
    output: "build/data/heat/annual-{building_sector}-{demand_type}.nc"
    script: "../scripts/heat/annual_heat_demand.py"

