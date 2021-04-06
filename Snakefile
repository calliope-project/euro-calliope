URL_LOAD = "https://data.open-power-system-data.org/time_series/2019-06-05/time_series_60min_stacked.csv"
URL_POTENTIALS = "https://zenodo.org/record/3533038/files/possibility-for-electricity-autarky.zip"
URL_CAPACITY_FACTORS = "https://zenodo.org/record/3899687/files/{wildcards.filename}?download=1"

NATIONAL_PHS_STORAGE_CAPACITIES = "data/pumped-hydro/storage-capacities-gwh.csv"
ALL_WIND_AND_SOLAR_TECHNOLOGIES = [
    "wind-onshore", "wind-offshore", "open-field-pv",
    "rooftop-pv", "rooftop-pv-n", "rooftop-pv-e-w", "rooftop-pv-s-flat"
]

include: "./rules/shapes.smk"
include: "./rules/hydro.smk"
localrules: all, raw_load, model, clean, parameterise_template, potentials_zipped
localrules: download_capacity_factors_wind_and_solar
configfile: "config/default.yaml"
wildcard_constraints:
        resolution = "((continental)|(national)|(regional))"

__version__ = open('./VERSION').readlines()[0].strip()

onstart:
    shell("mkdir -p build/logs")


rule all:
    message: "Generate euro-calliope pre-built models and run tests."
    input:
        "build/logs/continental/model.done",
        "build/logs/national/model.done",
        "build/logs/regional/model.done",
        "build/logs/continental/test-report.html",
        "build/logs/national/test-report.html",


rule potentials_zipped:
    message: "Download potential data."
    output: protected("data/automatic/raw-potentials.zip")
    shell: "curl -sLo {output} '{URL_POTENTIALS}'"


rule potentials:
    message: "Unzip potentials."
    input: rules.potentials_zipped.output[0]
    shadow: "minimal"
    output:
        land_eligibility_km2 = "build/data/{resolution}/technical-potential/areas.csv",
        shared_coast = "build/data/{resolution}/shared-coast.csv",
        industrial_demand = "build/data/{resolution}/demand.csv",
        population = "build/data/{resolution}/population.csv",
        land_cover = "build/data/{resolution}/land-cover.csv"
    shell: "unzip -o {input} -d build/data"


rule parameterise_template:
    message: "Apply config parameters to file {wildcards.template} from templates."
    input:
        src = "src/parameterise_templates.py",
        filters = "src/filters.py",
        template = "src/template/{template}",
        biofuel_cost = "build/data/regional/biofuel/{scenario}/costs-eur-per-mwh.csv".format(
            scenario=config["parameters"]["jrc-biofuel"]["scenario"]
        )
    params:
        scaling_factors = config["scaling-factors"],
        capacity_factors = config["capacity-factors"]["average"],
        max_power_density = config["parameters"]["maximum-installable-power-density"],
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"]
    output: "build/model/{template}"
    wildcard_constraints:
        template = "((link-techs.yaml)|(storage-techs.yaml)|(demand-techs.yaml)|(renewable-techs.yaml)|(README.md)|(environment.yaml)|(interest-rate.yaml))"
    conda: "envs/default.yaml"
    script: "src/parameterise_templates.py"


rule hydro_capacities:
    message: "Determine hydro capacities on {wildcards.resolution} resolution."
    input:
        src = "src/hydro.py",
        locations = rules.units.output[0],
        plants = rules.preprocess_hydro_stations.output[0],
        phs_storage_capacities = NATIONAL_PHS_STORAGE_CAPACITIES
    output: "build/data/{resolution}/hydro-capacities-mw.csv"
    conda: "envs/geo.yaml"
    script: "src/hydro.py"


BIOFUEL_FEEDSTOCKS = [
    "forestry-energy-residues",
    "landscape-care-residues",
    "manure",
    "municipal-waste",
    "primary-agricultural-residues",
    "roundwood-chips",
    "roundwood-fuelwood",
    "secondary-forestry-residues-sawdust",
    "secondary-forestry-residues-woodchips",
    "sludge"
]

rule biofuels:
    message: "Determine biofuels potential on {wildcards.resolution} resolution for scenario {wildcards.scenario}."
    input:
        src = "src/biofuels.py",
        units = rules.units.output[0],
        land_cover = rules.potentials.output.land_cover,
        population = rules.potentials.output.population,
        national_potentials = expand("data/biofuels/potentials/{feedstock}.csv", feedstock=BIOFUEL_FEEDSTOCKS),
        costs = expand("data/biofuels/costs/{feedstock}.csv", feedstock=BIOFUEL_FEEDSTOCKS)
    params:
        potential_year = config["parameters"]["jrc-biofuel"]["potential-year"],
        cost_year = config["parameters"]["jrc-biofuel"]["cost-year"]
    output:
        potentials = "build/data/{resolution}/biofuel/{scenario}/potential-mwh-per-year.csv",
        costs = "build/data/{resolution}/biofuel/{scenario}/costs-eur-per-mwh.csv" # not actually resolution dependent
    conda: "envs/geo.yaml"
    wildcard_constraints:
        scenario = "((low)|(medium)|(high))"
    script: "src/biofuels.py"


rule locations:
    message: "Generate locations for {wildcards.resolution} resolution."
    input:
        src = "src/locations.py",
        filters = "src/filters.py",
        shapes = rules.units.output[0],
        land_eligibility_km2 = rules.potentials.output.land_eligibility_km2,
        hydro_capacities = rules.hydro_capacities.output[0],
        biofuel = "build/data/{{resolution}}/biofuel/{scenario}/potential-mwh-per-year.csv".format(scenario=config["parameters"]["jrc-biofuel"]["scenario"])
    params:
        flat_roof_share = config["parameters"]["roof-share"]["flat"],
        maximum_installable_power_density = config["parameters"]["maximum-installable-power-density"],
        scaling_factors = config["scaling-factors"],
        biofuel_efficiency = config["parameters"]["biofuel-efficiency"]
    output:
        yaml = "build/model/{resolution}/locations.yaml",
        csv = "build/model/{resolution}/locations.csv"
    conda: "envs/geo.yaml"
    script: "src/locations.py"


rule directional_rooftop_pv:
    message: "Generate override for directional rooftop PV in {wildcards.resolution} resolution."
    input:
        src = "src/directional_rooftop.py",
        filters = "src/filters.py",
        shapes = rules.units.output[0],
        land_eligibility_km2 = rules.potentials.output.land_eligibility_km2,
    params:
        roof_shares = config["parameters"]["roof-share"],
        maximum_installable_power_density = config["parameters"]["maximum-installable-power-density"],
        scaling_factors = config["scaling-factors"],
    output: "build/model/{resolution}/directional-rooftop.yaml"
    conda: "envs/geo.yaml"
    script: "src/directional_rooftop.py"


rule load_shedding:
    message: "Generate override allowing load shedding."
    input:
        src = "src/load_shedding.py",
        shapes = rules.units.output[0]
    output: "build/model/{resolution}/load-shedding.yaml"
    conda: "envs/geo.yaml"
    script: "src/load_shedding.py"


rule download_capacity_factors_wind_and_solar:
    message: "Download data/automatic/capacityfactors/{wildcards.filename}."
    output: protected("data/automatic/capacityfactors/{filename}")
    shell: f"curl -sLo {{output}} '{URL_CAPACITY_FACTORS}'"


rule capacity_factors_onshore_wind_and_solar:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for {wildcards.technology}."
    input:
        src = "src/capacityfactors.py",
        locations = rules.units.output[0],
        ids = ancient("data/automatic/capacityfactors/onshore-locations.tif"),
        timeseries = ancient("data/automatic/capacityfactors/{technology}-timeseries.nc")
    params:
        threshold = config["capacity-factors"]["min"],
        year = config["year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    wildcard_constraints:
        technology = "((wind-onshore)|(rooftop-pv)|(open-field-pv)|(rooftop-pv-n)|(rooftop-pv-e-w)|(rooftop-pv-s-flat))"
    output: "build/model/{resolution}/capacityfactors-{technology}.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors.py"


rule capacity_factors_offshore:
    message: "Generate capacityfactor time series disaggregated by location on "
             "{wildcards.resolution} resolution for wind-offshore."
    input:
        src = "src/capacityfactors_offshore.py",
        eez = rules.eez.output[0],
        shared_coast = rules.potentials.output.shared_coast,
        ids = ancient("data/automatic/capacityfactors/offshore-locations.tif"),
        timeseries = ancient("data/automatic/capacityfactors/wind-offshore-timeseries.nc")
    params:
        threshold = config["capacity-factors"]["min"],
        year = config["year"],
        trim_ts = config["capacity-factors"]["trim-ninja-timeseries"]
    output: "build/model/{resolution}/capacityfactors-wind-offshore.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors_offshore.py"


rule capacity_factors_hydro:
    message: "Generate capacityfactor time series for hydro electricity on {wildcards.resolution} resolution."
    input:
        src = "src/capacityfactors_hydro.py",
        capacities = rules.hydro_capacities.output[0],
        stations = rules.inflow_mwh.output[0],
        locations = rules.units.output[0]
    params:
        threshold = config["capacity-factors"]["min"]
    output:
        ror = "build/model/{resolution}/capacityfactors-hydro-ror.csv",
        reservoir = "build/model/{resolution}/capacityfactors-hydro-reservoir-inflow.csv"
    conda: "envs/geo.yaml"
    script: "src/capacityfactors_hydro.py"


rule raw_load:
    message: "Download raw load."
    output: protected("data/automatic/raw-load-data.csv")
    shell: "curl -sLo {output} '{URL_LOAD}'"


rule electricity_load_national:
    message: "Preprocess raw electricity load data and retrieve load time series per country."
    input:
        src = "src/national_load.py",
        load = rules.raw_load.output
    output: "build/data/electricity-demand-national.csv"
    params:
        number_rows_valid = 20950674, # see https://github.com/Open-Power-System-Data/time_series/issues/22
        year = config["year"]
    conda: "envs/default.yaml"
    script: "src/national_load.py"


rule electricity_load:
    message: "Generate electricity load time series for every location on {wildcards.resolution} resolution."
    input:
        src = "src/load.py",
        units = rules.units.output[0],
        industrial_demand = rules.potentials.output.industrial_demand,
        national_load = rules.electricity_load_national.output[0]
    params:
        scaling_factor = config["scaling-factors"]["power"]
    output: "build/model/{resolution}/electricity-demand.csv"
    conda: "envs/geo.yaml"
    script: "src/load.py"


rule link_neighbours:
    message: "Create links between all direct neighbours on {wildcards.resolution} resolution."
    input:
        src = "src/link_neighbours.py",
        units = rules.units.output[0]
    params:
        sea_connections = lambda wildcards: config["sea-connections"][wildcards.resolution]
    output: "build/model/{resolution}/link-all-neighbours.yaml"
    conda: "envs/geo.yaml"
    script: "src/link_neighbours.py"


rule build_metadata:
    message: "Generate build metadata."
    input:
        src = "src/metadata.py"
    params:
        config = config,
        version = __version__
    output: "build/model/build-metadata.yaml"
    conda: "envs/default.yaml"
    script: "src/metadata.py"


rule model:
    message: "Generate euro-calliope with {wildcards.resolution} resolution."
    input:
        "build/model/interest-rate.yaml",
        "build/model/link-techs.yaml",
        "build/model/renewable-techs.yaml",
        "build/model/storage-techs.yaml",
        "build/model/demand-techs.yaml",
        "build/model/environment.yaml",
        "build/model/README.md",
        rules.locations.output,
        rules.electricity_load.output,
        rules.link_neighbours.output,
        rules.capacity_factors_hydro.output,
        rules.hydro_capacities.output,
        rules.directional_rooftop_pv.output,
        expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=ALL_WIND_AND_SOLAR_TECHNOLOGIES
        ),
        rules.build_metadata.output,
        example_model = "src/template/example-model.yaml"
    output:
        log = "build/logs/{resolution}/model.done",
        example_model = "build/model/{resolution}/example-model.yaml"
    shell:
        """
        cp {input.example_model} {output.example_model}
        touch {output.log}
        """


rule clean: # removes all generated results
    shell:
        """
        rm -r build/
        echo "Data downloaded to data/ has not been cleaned."
        """


rule test:
    message: "Run tests"
    input:
        "tests/test_runner.py",
        "tests/test_model.py",
        "tests/test_capacityfactors.py",
        "build/logs/{resolution}/model.done",
        model = "tests/resources/{resolution}/model.yaml",
        example_model = "build/model/{resolution}/example-model.yaml",
        capacity_factor_timeseries = expand(
            "build/model/{{resolution}}/capacityfactors-{technology}.csv",
            technology=ALL_WIND_AND_SOLAR_TECHNOLOGIES + ["hydro-ror", "hydro-reservoir-inflow"]
        )
    params:
        config = config
    output: "build/logs/{resolution}/test-report.html"
    conda: "./envs/test.yaml"
    script: "./tests/test_runner.py"
