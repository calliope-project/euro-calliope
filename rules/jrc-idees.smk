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

def get_countries_to_unzip(config: dict) -> set:
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
    message: "Unzip all JRC-IDEES {wildcards.sector} sector country data"
    input:
        countries = expand(
            "data/automatic/jrc-idees/{country_code}.zip",
            country_code=get_countries_to_unzip(config)
        )
    params:
        sector_title_case = lambda wildcards: wildcards.sector.title(),
        # These params are referenced so that the rule re-runs when they change value.
        # They are used in the `get_countries_to_unzip` helper function.
        infill_country_config = config["data-pre-processing"]["fill-missing-values"]["jrc-idees"],
        countries = config["scope"]["spatial"]["countries"]
    wildcard_constraints:
        sector = "industry|transport|tertiary"
    output: temp(directory("build/data/jrc-idees/{sector}/unprocessed"))
    conda: "../envs/shell.yaml"
    shell:
        """
        for file in {input.countries};
            do unzip "$file" '*{params.sector_title_case}*' -d {output};
        done
        """


rule jrc_idees_industry_processed:
    message: "Process {wildcards.dataset} industry data from JRC-IDEES to be used in understanding current and future industry demand"
    input:
        unprocessed_data = "build/data/jrc-idees/industry/unprocessed"
    output: "build/data/jrc-idees/industry/processed-{dataset}.nc"
    wildcard_constraints:
        dataset = "energy|production"
    conda: "../envs/default.yaml"
    threads: 4
    script: "../scripts/jrc-idees/industry.py"
