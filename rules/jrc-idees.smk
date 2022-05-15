"""Rules to process input data."""

localrules: jrc_idees_zipped

EU28 = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]


rule jrc_idees_zipped:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: url = config["data-sources"]["jrc-idees"]
    output: protected("data/automatic/jrc-idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLo {output} '{params.url}'"


rule jrc_idees_unzipped:
    message: "Unzip all JRC-IDEES {wildcards.sector} sector country data"
    input:
        countries = [
            f"data/automatic/jrc-idees/{country_code}.zip"
            for country_code in [
                pycountry.countries.lookup(country).alpha_2 for country in config["scope"]["spatial"]["countries"]
            ]
            if country_code in EU28
        ]
    params: sector_title_case = lambda wildcards: wildcards.sector.title()
    wildcard_constraints:
        sector = "((industry)|(transport)|(tertiary))"
    output: temp(directory("build/data/jrc-idees/{sector}/unprocessed"))
    conda: "../envs/shell.yaml"
    shell: "unzip 'data/automatic/jrc-idees/*.zip' '*{params.sector_title_case}*' -d {output}"


rule jrc_idees_transport_processed:
    message: "Process {wildcards.dataset} transport data from JRC-IDEES to be used in understanding current and future transport demand"
    input:
        script = script_dir + "jrc-idees/transport.py",
        unprocessed_data = "build/data/jrc-idees/transport/unprocessed"
    output: "build/data/jrc-idees/transport/processed-{dataset}.csv"
    wildcard_constraints:
        dataset = "((road-energy)|(road-distance)|(road-vehicles)|(rail-energy)|(rail-distance))"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/transport.py"


rule jrc_idees_industry_processed:
    message: "Process {wildcards.dataset} industry data from JRC-IDEES to be used in understanding current and future industry demand"
    input:
        script = script_dir + "jrc-idees/industry.py",
        unprocessed_data = "build/data/jrc-idees/industry/unprocessed"
    output: "build/data/jrc-idees/industry/processed-{dataset}.csv"
    wildcard_constraints:
        dataset = "((energy)|(production))"
    conda: "../envs/default.yaml"
    threads: config["snakemake"]["max-threads"]
    script: "../scripts/jrc-idees/industry.py"


rule jrc_idees_tertiary_processed:
    message: "Process tertiary sector energy data from JRC-IDEES to be used in understanding current and future tertiary sector demand"
    input:
        script = script_dir + "jrc-idees/tertiary.py",
        unprocessed_data = "build/data/jrc-idees/tertiary/unprocessed"
    output: "build/data/jrc-idees/tertiary/processed-energy.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/tertiary.py"
