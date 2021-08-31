"""Rules that access raw data."""

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
localrules: download_euro_calliope_datasets, euro_calliope_datasets


rule download_euro_calliope_datasets:
    message: "Download Euro-Calliope datasets."
    params: url = config["data-sources"]["data-repository"]
    output: protected("data/automatic/euro-calliope-datasets.zip")
    conda: "../envs/shell.yaml"
    shell: "curl -sLfo {output} {params.url}"


rule euro_calliope_datasets:
    message: "Unzip Euro-Calliope datasets."
    input: rules.download_euro_calliope_datasets.output[0]
    shadow: "minimal"
    output:
        irena_generation = "build/data/irena/hydro-generation-europe.csv",
        national_phs_storage_capacities = "build/data/pumped-hydro/storage-capacities-gwh.csv",
        biofuel_potentials = expand(
            "build/data/biofuels/potentials/{feedstock}.csv",
            feedstock=BIOFUEL_FEEDSTOCKS),
        biofuel_costs = expand(
            "build/data/biofuels/costs/{feedstock}.csv",
            feedstock=BIOFUEL_FEEDSTOCKS)
    conda: "../envs/shell.yaml"
    shell:
        # This script removes the top-level directory, whose name we do not know.
        # Because tar can do this but unzip can not, we unzip first, then tar again,
        # and finally untar while removing the top-level.
        """
        mkdir build/tmp/
        unzip -o {input} -d build/tmp
        tar cf build/tmp.tar --directory build/tmp/ .
        tar xf build/tmp.tar --directory build/data --strip-components=2
        """
