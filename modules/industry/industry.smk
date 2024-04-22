# Paths dependent on main Snakefile
MODULE_PATH = "modules/industry"
BUILD_PATH = f"{MODULE_PATH}/build"
DATA_PATH = f"{MODULE_PATH}/raw_data"

# Paths relative to this snakefile (snakemake behaviour is inconsitent)
SCRIPT_PATH = "scripts"  # scripts are called relative to this file
CONDA_PATH = "./env_industry.yaml"

# TODO:
# jrc_idees_processed* files in the raw_data folder should ideally be produced by a rule instead

# Ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work!
rule steel_industry:
    message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
    conda: CONDA_PATH
    params:
        cnf_steel = config["params"]["steel"]
    input:
        path_energy_balances = config["inputs"]["path-energy-balances"],
        path_cat_names = config["inputs"]["path-cat-names"],
        path_carrier_names = config["inputs"]["path-carrier-names"],
        path_jrc_industry_energy = config["inputs"]["path-jrc-industry-energy"],
        path_jrc_industry_production = config["inputs"]["path-jrc-industry-production"],
    output:
        path_output = f"{BUILD_PATH}/annual_demand_steel.nc"
    script: f"{SCRIPT_PATH}/steel_industry.py"

rule chemical_industry:
    message: "."
    conda: CONDA_PATH
    params:
    input:
    output:
    script: f"{SCRIPT_PATH}/chemicals.py"

rule other_industry:
    message: "."
    conda: CONDA_PATH
    params:
    input:
    output: f"{BUILD_PATH}/other_industry.csv"
    script: f"{SCRIPT_PATH}/other_industry.py"

# rule combine_and_scale:
#     message: "."
#     conda: CONDA_PATH
#     params:
#     input:
#     output:
#     script:

# rule verify:
#     message: "."
#     params:
#     input:
#     output:
#     script:
