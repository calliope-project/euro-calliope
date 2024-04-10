SRC_PATH = "src"
TMP_PATH = "modules/industry/tmp"
OUT_PATH = "modules/industry/out"
DATA_PATH = "modules/industry/raw_data"
# Snakemake searches for conda environments at the same level as this file
CONDA_PATH = "./env_industry.yaml"

# TODO:
# jrc_idees_processed* files in the raw_data folder should ideally be produced by a rule instead

# Ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work!
rule steel_industry:
    message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
    conda: CONDA_PATH
    params:
        year_range = config["params"]["year-range"],
        path_energy_balances = config["inputs"]["path-energy-balances"],
        path_cat_names = config["inputs"]["path-cat-names"],
        path_carrier_names = config["inputs"]["path-carrier-names"],
    input:
        path_jrc_energy = f"{DATA_PATH}/jrc_idees_processed_energy.csv.gz",
        path_jrc_production = f"{DATA_PATH}/jrc_idees_processed_production.csv.gz",
    output:
        path_output = f"{TMP_PATH}/annual_demand_steel.csv"
    script: f"{SRC_PATH}/steel_industry.py"

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
    message: "."
    conda: CONDA_PATH
    params:
    input:
    output: f"{TMP_PATH}/other_industry.csv"
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
