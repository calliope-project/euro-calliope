from snakemake.utils import validate

# Paths dependent on main Snakefile
MODULE_PATH = "modules/industry"  # TODO: remove if the module becomes an imported external workflow
BUILD_PATH = f"{MODULE_PATH}/build"
DATA_PATH = f"{MODULE_PATH}/raw_data"

# Paths relative to this snakefile (snakemake behaviour is inconsitent)
SCRIPT_PATH = "scripts"  # scripts are called relative to this file
UTILS = [f"{MODULE_PATH}/{SCRIPT_PATH}/utils/{i}.py" for i in ["filling", "jrc_idees_parser"]]
CONDA_PATH = "./env_industry.yaml"

configfile: "./config.yaml"
validate(config, "./schema.yaml")

# Ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work!
rule iron_and_steel:
    message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
    conda: CONDA_PATH
    params:
        config = config["params"]["config-iron-and-steel"]
    input:
        energy_balances = config["input-paths"]["energy-balances"],
        cat_names = config["input-paths"]["cat-names"],
        carrier_names = config["input-paths"]["carrier-names"],
        jrc_industry_energy = config["input-paths"]["jrc-industry-energy"],
        jrc_industry_production = config["input-paths"]["jrc-industry-production"],
    output:
        path_output = f"{BUILD_PATH}/annual_demand_iron_and_steel.nc"
    script: f"{SCRIPT_PATH}/iron_and_steel.py"

rule chemicals_industry:
    message: "Calculate energy demand for the 'Chemicals Industry' sector in JRC-IDEES."
    conda: CONDA_PATH
    # params:
    #     config = config["params"]["config-chemicals-industry"]
    input:
        energy_balances = config["input-paths"]["energy-balances"],
        cat_names = config["input-paths"]["cat-names"],
        carrier_names = config["input-paths"]["carrier-names"],
        jrc_industry_energy = config["input-paths"]["jrc-industry-energy"],
        jrc_industry_production = config["input-paths"]["jrc-industry-production"],
    output:
        path_output = f"{BUILD_PATH}/annual_demand_chemicals_industry.nc"
    script: f"{SCRIPT_PATH}/chemicals_industry.py"

rule combined_categories:
    message: "Calculate energy demand for all other industry sectors in JRC-IDEES."
    conda: CONDA_PATH
    params:
        specific_categories = config["params"]["specific-categories"],
        config = config["params"]["config-combined-categories"],
    input:
        energy_balances = config["input-paths"]["energy-balances"],
        cat_names = config["input-paths"]["cat-names"],
        carrier_names = config["input-paths"]["carrier-names"],
        jrc_industry_energy = config["input-paths"]["jrc-industry-energy"],
        jrc_industry_production = config["input-paths"]["jrc-industry-production"],
    output: f"{BUILD_PATH}/annual_demand_combined_categories.nc"
    script: f"{SCRIPT_PATH}/combined_categories.py"

SUFFIXES = [i.lower().replace(" ", "_") for i in config["params"]["specific-categories"]]
rule combine_and_scale:
    message: "Identify the category scripts to run based on the configuration."
    conda: CONDA_PATH
    input:
        UTILS,
        expand("{path}/annual_demand_{sample}.nc", path=[BUILD_PATH], sample=SUFFIXES),
        rules.combined_categories.output
    output: "{BUILD_PATH}/annual_demand_aggregated.nc"
    shell: "touch {output}"


# rule verify:
#     message: "."
#     params:
#     input:
#     output:
#     script:
