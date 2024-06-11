from snakemake.utils import validate

# Paths dependent on main Snakefile
MODULE_PATH = "modules/industry"
BUILD_PATH = f"{MODULE_PATH}/build"
DATA_PATH = f"{MODULE_PATH}/raw_data"

# Paths relative to this snakefile (snakemake behaviour is inconsitent)
SCRIPT_PATH = "scripts"  # scripts are called relative to this file
CONDA_PATH = "./env_industry.yaml"

validate(config, "./schema.yaml")

# Ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work!
if "Iron and steel" in config["params"]["non-generic-categories"]:
    rule steel_processing:
        message: "Calculate energy demand for the 'Iron and steel' sector in JRC-IDEES."
        conda: CONDA_PATH
        params:
            steel_config = config["params"]["steel-config"]
        input:
            path_energy_balances = config["inputs"]["path-energy-balances"],
            path_cat_names = config["inputs"]["path-cat-names"],
            path_carrier_names = config["inputs"]["path-carrier-names"],
            path_jrc_industry_energy = config["inputs"]["path-jrc-industry-energy"],
            path_jrc_industry_production = config["inputs"]["path-jrc-industry-production"],
        output:
            path_output = f"{BUILD_PATH}/annual_demand_steel.nc"
        script: f"{SCRIPT_PATH}/steel_processing.py"

if "Chemicals Industry" in config["params"]["non-generic-categories"]:
    rule chemicals_processing:
        message: "."
        conda: CONDA_PATH
        params:
        input:
        output:
        script: f"{SCRIPT_PATH}/chemicals_processing.py"

rule generic_processing:
    message: "Calculate energy demand for all other industry sectors in JRC-IDEES."
    conda: CONDA_PATH
    params:
        non_generic_categories = config["params"]["non-generic-categories"],
        generic_config = config["params"]["generic-config"],
    input:
        path_energy_balances = config["inputs"]["path-energy-balances"],
        path_cat_names = config["inputs"]["path-cat-names"],
        path_carrier_names = config["inputs"]["path-carrier-names"],
        path_jrc_industry_energy = config["inputs"]["path-jrc-industry-energy"],
        path_jrc_industry_production = config["inputs"]["path-jrc-industry-production"],
    output:
        path_output = f"{BUILD_PATH}/annual_demand_generic.nc"
    script: f"{SCRIPT_PATH}/generic_processing.py"

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
