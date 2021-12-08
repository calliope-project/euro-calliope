"""Rules related to biofuels."""

root_dir = config["root-directory"] + "/" if config["root-directory"] not in ["", "."] else ""
script_dir = f"{root_dir}scripts/"


rule download_biofuel_potentials_and_costs:
    message: "Download raw biofuel potential and cost data."
    params: url = config["data-sources"]["biofuel-potentials-and-costs"]
    output: protected("data/automatic/raw-biofuel-potentials-and-costs.xslx")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule preprocess_biofuel_potentials_and_cost:
    message: "Extract national potentials and cost from raw biofuel data."
    input:
        script = script_dir + "preprocess/extract_biofuels.py",
        potentials_and_costs = rules.download_biofuel_potentials_and_costs.output[0]
    output:
        potentials = "build/data/raw-biofuel-potentials.csv",
        costs = "build/data/raw-biofuel-costs.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/preprocess/extract_biofuels.py"


rule biofuels:
    message: "Determine biofuels potential on {wildcards.resolution} resolution for scenario {wildcards.scenario}."
    input:
        script = script_dir + "biofuels.py",
        units = rules.units_without_shape.output[0],
        land_cover = rules.potentials.output.land_cover,
        population = rules.potentials.output.population,
        national_potentials = rules.preprocess_biofuel_potentials_and_cost.output.potentials,
        costs = rules.preprocess_biofuel_potentials_and_cost.output.costs
    params:
        potential_year = config["parameters"]["jrc-biofuel"]["potential-year"],
        cost_year = config["parameters"]["jrc-biofuel"]["cost-year"]
    output:
        potentials = "build/data/{resolution}/biofuel/{scenario}/potential-mwh-per-year.csv",
        costs = "build/data/{resolution}/biofuel/{scenario}/costs-eur-per-mwh.csv" # not actually resolution dependent
    conda: "../envs/default.yaml"
    wildcard_constraints:
        scenario = "((low)|(medium)|(high))"
    script: "../scripts/biofuels.py"
