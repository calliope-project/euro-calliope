import pandas as pd


def biomass_potential_mw(path_to_population, path_to_farmland, annual_biomass_from_waste_per_capita,
                         biomass_yield_from_crops, path_to_output):
    """Determine the energy stream of biomass available in each location."""
    hourly_waste_energy = potential_from_waste(path_to_population, annual_biomass_from_waste_per_capita)
    hourly_crop_energy = potential_from_crops(path_to_farmland, biomass_yield_from_crops)
    hourly_energy = hourly_waste_energy + hourly_crop_energy
    hourly_energy.name = "biomass_potential_mwh_per_hour"
    validate_result(hourly_energy)
    hourly_energy.to_csv(
        path_to_output,
        index=True,
        header=True
    )


def potential_from_waste(path_to_population, annual_biomass_from_waste_per_capita):
    population_count = pd.read_csv(path_to_population, index_col=0)["population_sum"]
    annual_waste_energy = population_count * annual_biomass_from_waste_per_capita
    return annual_waste_energy / 365 / 24


def potential_from_crops(path_to_farmland, biomass_yield_from_crops):
    # ASSUME wind and biomass can be used together without restrictions
    # ASSUME not eligble farmland is not eligible for wind and pv because of the slope (true for technical-potential)
    # ASSUME all farmland can be used for bioenergy
    # ASSUME pv and biomass can be used together without restrictions
    # TODO ^ revisit the last two assumptions which are problematic
    farmland = pd.read_csv(path_to_farmland, index_col=0)
    return (
        farmland["eligibility_not_eligible_km2"] * biomass_yield_from_crops
        + farmland["eligibility_onshore_wind_km2"] * biomass_yield_from_crops
        + farmland["eligibility_onshore_wind_and_pv_km2"] * biomass_yield_from_crops
    )


def validate_result(hourly_energy):
    # technical potential of biomass is as high as 50% to 100% of current electricity demand
    n_hours_per_year = 8760
    european_demand = 3_500_000_000 # [MWh], roughly
    conversion_efficiency = 0.4 # roughly
    annual_usable_energy = hourly_energy * n_hours_per_year * conversion_efficiency
    assert european_demand / 2 < annual_usable_energy.sum() < european_demand


if __name__ == "__main__":
    biomass_potential_mw(
        path_to_population=snakemake.input.population,
        path_to_farmland=snakemake.input.farmland,
        annual_biomass_from_waste_per_capita=snakemake.params.annual_biomass_from_waste_per_capita,
        biomass_yield_from_crops=snakemake.params.biomass_yield_from_crops,
        path_to_output=snakemake.output[0]
    )
