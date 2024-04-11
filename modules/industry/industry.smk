# Paths dependent on main Snakefile
MODULE_PATH = "modules/industry"
TMP_PATH = f"{MODULE_PATH}/tmp"
OUT_PATH = f"{MODULE_PATH}/out"
DATA_PATH = f"{MODULE_PATH}/raw_data"

# Paths relative to this snakefile (snakemake behaviour is inconsitent)
SRC_PATH = "src"  # scripts are called relative to this file
CONDA_PATH = "./env_industry.yaml"

# TODO:
# jrc_idees_processed* files in the raw_data folder should ideally be produced by a rule instead

# Ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work!
if "Iron and steel" in config["params"]["specific-industries"]:
    rule steel_industry:
        message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
        conda: CONDA_PATH
        params:
            year_range = config["params"]["year-range"]
        input:
            path_energy_balances = config["inputs"]["path-energy-balances"],
            path_cat_names = config["inputs"]["path-cat-names"],
            path_carrier_names = config["inputs"]["path-carrier-names"],
            path_jrc_energy = f"{DATA_PATH}/jrc_idees_processed_energy.csv.gz",
            path_jrc_production = f"{DATA_PATH}/jrc_idees_processed_production.csv.gz",
        output:
            path_output = f"{TMP_PATH}/annual_demand_steel.csv"
        script: f"{SRC_PATH}/steel_industry.py"

if "Chemicals Industry" in config["params"]["specific-industries"]:
    rule chemicals_industry:
        message: "Calculate energy demand for the 'Chemicals Industry' sector in JRC-IDEES."
        conda: CONDA_PATH
        params:
            year_range = config["params"]["year-range"],
        input:
            path_energy_balances = config["inputs"]["path-energy-balances"],
            path_cat_names = config["inputs"]["path-cat-names"],
            path_carrier_names = config["inputs"]["path-carrier-names"],
            path_jrc_energy = f"{DATA_PATH}/jrc_idees_processed_energy.csv.gz",
            path_jrc_production = f"{DATA_PATH}/jrc_idees_processed_production.csv.gz",
        output:
            path_output = f"{TMP_PATH}/annual_demand_chemicals.csv"
        script: f"{SRC_PATH}/chemicals_industry.py"

rule other_industry:
    message: "Calculate energy demand for all other industry sectors in JRC-IDEES."
    conda: CONDA_PATH
    params:
        year_range = config["params"]["year-range"],
        specific_industries = config["params"]["specific-industries"]
    input:
        path_energy_balances = config["inputs"]["path-energy-balances"],
        path_cat_names = config["inputs"]["path-cat-names"],
        path_carrier_names = config["inputs"]["path-carrier-names"],
        path_jrc_energy = f"{DATA_PATH}/jrc_idees_processed_energy.csv.gz",
        path_jrc_production = f"{DATA_PATH}/jrc_idees_processed_production.csv.gz",
    output:
        path_output = f"{TMP_PATH}/annual_demand_other.csv"
    script: f"{SRC_PATH}/other_industry.py"

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
