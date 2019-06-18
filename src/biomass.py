import pandas as pd


def biomass_potential_mw(path_to_population, annual_biomass_from_waste_per_capita, path_to_output):
    """Determine the energy stream of biomass available in each location."""
    population_count = pd.read_csv(path_to_population, index_col=0)["population_sum"]
    annual_waste_energy = population_count * annual_biomass_from_waste_per_capita
    hourly_waste_energy = annual_waste_energy / 365 / 24
    hourly_waste_energy.name = "biomass_potential_mwh_per_hour"
    hourly_waste_energy.to_csv(
        path_to_output,
        index=True,
        header=True
    )


if __name__ == "__main__":
    biomass_potential_mw(
        path_to_population=snakemake.input.population,
        annual_biomass_from_waste_per_capita=snakemake.params.annual_biomass_from_waste_per_capita,
        path_to_output=snakemake.output[0]
    )
