"Rules regarding JRC-IDEES Data"

JRC_IDEES_SPATIAL_SCOPE = [
    "AT", "BE", "BG", "CY", "CZ", "DE", "DK", "EE", "EL", "ES", "FI", "FR",
    "HR", "HU", "IE", "IT", "LT", "LU", "LV", "MT", "NL", "PL", "PT", "RO",
    "SE", "SI", "SK", "UK"
]


rule download_jrc_idees_zipped:
    message: "Download JRC IDEES zip file for {wildcards.country_code}"
    params: url = config["data-sources"]["jrc-idees"]
    output: protected("data/automatic/jrc-idees/{country_code}.zip")
    conda: "../envs/shell.yaml"
    localrule: True
    shell: "curl -sSLo {output} '{params.url}'"


rule jrc_idees_unzipped:
    message: "Unzip JRC-IDEES {wildcards.sector} sector data for {wildcards.country_code}"
    input:
        country_data = "data/automatic/jrc-idees/{country_code}.zip",
    params:
        sector_title_case = lambda wildcards: wildcards.sector.title()
    wildcard_constraints:
        sector = "industry|transport|tertiary"
    output: temp("build/data/jrc-idees/{sector}/unprocessed/{country_code}.xlsx")
    conda: "../envs/shell.yaml"
    shell: "unzip -p {input.country_data} JRC-IDEES-2015_{params.sector_title_case}_{wildcards.country_code}.xlsx > {output}"


rule jrc_idees_industry_processed:
    message: "Process {wildcards.dataset} industry data from JRC-IDEES to be used in understanding current and future industry demand"
    input:
        data = expand(
            "build/data/jrc-idees/industry/unprocessed/{country_code}.xlsx",
            country_code=JRC_IDEES_SPATIAL_SCOPE
        )
    output: "build/data/jrc-idees/industry/processed-{dataset}.nc"
    wildcard_constraints:
        dataset = "energy|production"
    conda: "../envs/default.yaml"
    threads: 4
    script: "../scripts/jrc-idees/industry.py"


rule jrc_idees_tertiary_processed:
    message: "Process tertiary heat data from JRC-IDEES"
    input:
        data = expand(
            "build/data/jrc-idees/tertiary/unprocessed/{country_code}.xlsx",
            country_code=JRC_IDEES_SPATIAL_SCOPE
        )
    output: "build/data/jrc-idees/tertiary/processed.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/heat.py"


rule jrc_idees_transport_processed:
    message: "Process {wildcards.dataset} transport data from JRC-IDEES to be used in understanding current and future transport demand"
    input:
        data = expand(
            "build/data/jrc-idees/transport/unprocessed/{country_code}.xlsx",
            country_code=JRC_IDEES_SPATIAL_SCOPE
        )
    output: "build/data/jrc-idees/transport/processed-{dataset}.csv"
    params:
        vehicle_type_names = config["parameters"]["transport"]["vehicle-type-names"],
    wildcard_constraints:
        dataset = "road-energy|road-distance|road-vehicles"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/transport.py"
