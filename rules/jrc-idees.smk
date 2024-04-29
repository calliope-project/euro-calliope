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

def get_countries_to_unzip(countries: list[str], missing_value_mapping: dict[str, list]) -> set:
    """Get all countries whose JRC-IDEES files are to be unzipped.

    Includes primary configured countries and any neighbours of those countries required for data infilling.
    """
    primary_countries = [
        pycountry.countries.lookup(country).alpha_2
        for country in config["scope"]["spatial"]["countries"]
    ]

    infill_country_config = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
    infill_countries = [
        country for country in primary_countries
        for infill_country in infill_country_config.get(pycountry.countries.lookup(country).alpha_3, [])
    ]

    # Convert ISO standard country codes to their EU equivalents
    all_countries = [
        i.replace("GB", "UK").replace("GR", "EL")
        for i in set(primary_countries).union(infill_countries)
    ]

    return set(all_countries).intersection(JRC_IDEES_SPATIAL_SCOPE)

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
            country_code=get_countries_to_unzip(
                config["scope"]["spatial"]["countries"],
                config["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
            )
        )
    output: "build/data/jrc-idees/industry/processed-{dataset}.nc"
    wildcard_constraints:
        dataset = "energy|production"
    conda: "../envs/default.yaml"
    threads: 4
    script: "../scripts/jrc-idees/industry.py"


rule jrc_idees_heat_processed:
    message: "Process tertiary heat data from JRC-IDEES"
    input:
        data = expand(
            "build/data/jrc-idees/heat/unprocessed/{country_code}.xlsx",
            country_code=get_countries_to_unzip(
                config["scope"]["spatial"]["countries"],
                config["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
            )
        )
    output: "build/data/jrc-idees/heat/tertiary/processed.csv"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/heat.py"


rule jrc_idees_transport_processed:
    message: "Process {wildcards.dataset} transport data from JRC-IDEES to be used in understanding current and future transport demand"
    input:
        data = expand(
            "build/data/jrc-idees/transport/unprocessed/{country_code}.xlsx",
            country_code=get_countries_to_unzip(
                config["scope"]["spatial"]["countries"],
                config["data-pre-processing"]["fill-missing-values"]["jrc-idees"]
            )
        )
    output: "build/data/jrc-idees/transport/processed-{dataset}.csv"
    params:
        vehicle_type_names = config["parameters"]["transport"]["vehicle-type-names"],
    wildcard_constraints:
        dataset = "road-energy|road-distance|road-vehicles"
    conda: "../envs/default.yaml"
    script: "../scripts/jrc-idees/transport.py"
