SRC_PATH = "src"
TMP_PATH = "modules/industry/tmp"
OUT_PATH = "modules/industry/out"
DATA_PATH = "modules/industry/raw_data"
CONDA_PATH = "./env_industry.yaml"

# You must ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work.
# See https://github.com/snakemake/snakemake/issues/2514
rule steel_iron_industry:
    message: "."
    conda: CONDA_PATH
    params:
    input:
        path_energy_balances = config["inputs"]["path-energy-balances"],
        path_jrc_energy = f"{DATA_PATH}/jrc_idees_processed_energy.csv.gz",
        path_jrc_production = f"{DATA_PATH}/jrc_idees_processed_production.csv.gz"
    output: f"{TMP_PATH}/annual_iron_steel_demand.csv"
    script: f"{SRC_PATH}/steel_iron_industry.py"

rule chemical_industry:
    message: "."
    params:
    input:
    output: f"{TMP_PATH}/output.csv"
    script: f"{SRC_PATH}/chemicals.py"

rule other_industry:
    message: "."
    params:
    input:
    output: f"{TMP_PATH}/other_industry.csv"
    shell: f"{SRC_PATH}/other_industry.py"

rule combine:
    message: "."
    params:
    input:
    output: f"{TMP_PATH}/local.txt"
    script: "touch {output}"

rule test:
    message: "."
    params:
    input:
    output: f"{TMP_PATH}/local.txt"
    script: "touch {output}"

rule output:
    message: "."
    params:
    input:
        rules.chemical_industry.output,
        rules.steel_iron_industry.output,
        rules.other_industry.output
