"""Rules to generate electricity demand time series."""
import pycountry

localrules: download_raw_load, download_population_count, unzip_population_count, download_national_energy_balances, unzip_national_energy_balances


rule download_raw_load:
    message: "Download raw load."
    params: url = config["data-sources"]["load"]
    output: protected("data/automatic/raw-load-data.csv")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule download_population_count:
    message: "Download population count as raster."
    params: url = config["data-sources"]["population"]
    output: protected("data/automatic/raw-population-count.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_population_count:
    message: "Unzip population count."
    input: rules.download_population_count.output
    output: "build/data/population/GHS_POP_E2015_GLOBE_R2019A_4326_30ss_V1_0.tif"
    shadow: "minimal"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} -d build/data/population"


rule add_population_to_units:
    message: "Add population count to units of resolution {wildcards.resolution}."
    input:
        units = rules.units.output,
        population = rules.unzip_population_count.output
    output: "build/data/{resolution}/units_with_population.geojson"
    conda: "../envs/geo.yaml"
    shell:
        """
        fio cat {input.units} | \
        rio zonalstats -r {input.population} --prefix 'population_' --stats sum > \
        {output}
        """


rule add_population_share_to_units:
    message: "Add the unit's share of national population, for all units of resolution {wildcards.resolution}."
    input:
        units = rules.add_population_to_units.output[0]
    output: "build/data/{resolution}/units_with_population_share.geojson"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/add_population_share_to_unit.py"


rule download_national_energy_balances:
    message: "Download national energy balances which include demand by industrial sector."
    params: url = config["data-sources"]["national-energy-balances"]
    output: protected("data/automatic/raw-energy-balances.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule unzip_national_energy_balances:
    message: "Unzip national energy balances which include demand by industrial sector."
    input: rules.download_national_energy_balances.output
    output:
        national_energy_balances = [f"build/data/energy-balances/{country_code}-Energy-balance-sheets-June-2021-edition.xlsb"
                                    for country_code in ['EL' if code == 'GR' else code for code in
                                                        ['UK' if code == 'GB' else code for code in
                                                        [pycountry.countries.lookup(country).alpha_2 for country in
                                                        config["scope"]["spatial"]["countries"]]]]
                                   ] # Since alpha 2 is not correct for GB (UK) and GR (EL).
    shadow: "minimal"
    conda: "../envs/shell.yaml"
    shell: "unzip {input} -d build/data/energy-balances"


rule extract_ind_elec_demand:
    message: "Extracts national electricity demand by sector from .xslb-energy balances to .csv."
    input:
        national_energy_balances = rules.unzip_national_energy_balances.output
    params:
        path_energy_balances_foldername = "build/data/energy-balances/",
        path_energy_balances_filename = "-Energy-balance-sheets-June-2021-edition.xlsb",
        countries = config["scope"]["spatial"]["countries"],
        year = "2019" # ASSUME ratio emissions of sector j in unit i / national emissions of sector j constant over years
    output: "build/data/energy-balances/nat-ind-elec-demand.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/demand/nat-ind-elec-demand.py"


rule compute_emission_fractions:
    message: "Compute the unit's share of the national pollutant emissions / electricity demand of each sector."
    input:
        industrial_emission_data_master = "data/industrial-emission-data-master.csv",
        units = rules.add_population_share_to_units.output[0],
        mapping_unit_codes_nat_codes = rules.units_without_shape.output[0]
    output: "build/data/{resolution}/demand/units-share-of-nat-demand-per-sector.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/compute_emission_fractions.py"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        script = script_dir + "demand/national_load.py",
        load = rules.download_raw_load.output[0]
    params:
        first_year = config["scope"]["temporal"]["first-year"],
        final_year = config["scope"]["temporal"]["final-year"],
        data_quality_config = config["quality-control"]["load"],
        countries = config["scope"]["spatial"]["countries"]
    output: "build/data/electricity-demand-national.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/demand/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        script = script_dir + "demand/load.py",
        units = rules.units.output[0],
        demand_per_unit = rules.potentials.output.demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/timeseries/demand/electricity.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/load.py"


rule electricity_load_NEW:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        units_with_population_share = rules.add_population_share_to_units.output[0],
        units_share_of_nat_demand_per_sector = rules.compute_emission_fractions.output[0],
        nat_ind_elec_demand = rules.extract_ind_elec_demand.output[0],
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor_power = config["scaling-factors"]["power"]
    output: "build/models/{resolution}/timeseries/demand/electricity_NEW.csv"
    conda: "../envs/geo.yaml"
    script: "../scripts/demand/load_NEW.py"
