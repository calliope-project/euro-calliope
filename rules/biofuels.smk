"""Rules related to biofuels."""

localrules: download_biofuel_potentials_and_costs


rule download_biofuel_potentials_and_costs:
    message: "Download raw biofuel potential and cost data."
    params: url = config["data-sources"]["biofuel-potentials-and-costs"]
    output: protected("data/automatic/raw-biofuel-potentials-and-costs.xlsx")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule preprocess_biofuel_potentials_and_cost:
    message: "Extract national potentials and cost from raw biofuel data."
    input:
        potentials_and_costs = rules.download_biofuel_potentials_and_costs.output[0]
    params:
        feedstocks = {
            feedstock["id"]: name
            for name, feedstock in config["parameters"]["jrc-biofuel"]["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "build/data/raw-biofuel-potentials.csv",
        costs = "build/data/raw-biofuel-costs.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/biofuels/extract.py"


rule biofuels:
    message: "Determine biofuels potential on {wildcards.resolution} resolution for scenario {wildcards.scenario}."
    input:
        units = rules.units_without_shape.output[0],
        land_cover = rules.potentials.output.land_cover,
        population = rules.potentials.output.population,
        national_potentials = rules.preprocess_biofuel_potentials_and_cost.output.potentials,
        costs = rules.preprocess_biofuel_potentials_and_cost.output.costs
    params:
        potential_year = config["parameters"]["jrc-biofuel"]["potential-year"],
        cost_year = config["parameters"]["jrc-biofuel"]["cost-year"],
        proxies = {
            name: feedstock["proxy"]
            for name, feedstock in config["parameters"]["jrc-biofuel"]["feedstocks"].items()
            if feedstock["include"]
        }
    output:
        potentials = "build/data/{resolution}/biofuel/{scenario}/potential-mwh-per-year.csv",
        costs = "build/data/{resolution}/biofuel/{scenario}/costs-eur-per-mwh.csv" # not actually resolution dependent
    conda: "../envs/default.yaml"
    wildcard_constraints:
        scenario = "low|medium|high"
    script: "../scripts/biofuels/allocate.py"


rule bio_techs_and_locations_template:
    message: "Create biofuel tech definition file from template."
    input:
        template = techs_template_dir + "supply/biofuel.yaml",
        biofuel_cost = "build/data/regional/biofuel/{scenario}/costs-eur-per-mwh.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        ),
        locations = "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
    params:
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"],
        scaling_factors = config["scaling-factors"],
    conda: "../envs/default.yaml"
    output: "build/models/{resolution}/techs/supply/biofuel.yaml"
    script: "../scripts/biofuels/template_bio.py"
