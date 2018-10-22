PYTHON = "PYTHONPATH=./ python"

TECHNICAL_POTENTIAL_COUNTRIES = "src/data/national-technical-potential.geojson"

rule all:
    message: "Generate Euro Calliope."
    input:
        "model/countries.yaml"


rule countries:
    message: "Generate locations for all countries."
    input: TECHNICAL_POTENTIAL_COUNTRIES
    output: "model/countries.yaml"
    conda: "src/envs/geo.yaml"
    script: "src/locations.py"


rule clean: # removes all generated results
    shell:
        """
        rm -r model/
        """


rule test:
    shell:
        "py.test"
