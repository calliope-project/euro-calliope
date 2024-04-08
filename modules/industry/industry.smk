DATAPATH = config["setup"]["data-path"]
OUTPATH = config["setup"]["out-path"]
SCRIPTPATH = config["setup"]["script-path"]

# You must ensure rules are defined in order.
# Otherwise commands like "rules.rulename.output" won't work.
# See https://github.com/snakemake/snakemake/issues/2514
rule chemical_industry:
    message: "."
    params:
    input:
    output: f"{DATAPATH}/chemical_industry.csv"
    script: f"{SCRIPTPATH}/chemical_industry.py"

rule steel_iron_industry:
    message: "."
    params:
    input:
    output: f"{DATAPATH}/steel_iron_industry.csv"
    shell: f"{SCRIPTPATH}/steel_iron_industry.py"

rule other_industry:
    message: "."
    params:
    input:
    output: f"{DATAPATH}/other_industry.csv"
    shell: f"{SCRIPTPATH}/other_industry.py"

rule combine:
    message: "."
    params:
    input:
    output: f"{DATAPATH}/local.txt"
    script: "touch {output}"

rule test:
    message: "."
    params:
    input:
    output: f"{DATAPATH}/local.txt"
    script: "touch {output}"

rule output:
    message: "."
    params:
    input:
        rules.chemical_industry.output,
        rules.steel_iron_industry.output,
        rules.other_industry.output
