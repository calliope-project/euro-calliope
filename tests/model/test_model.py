import pandas as pd
import pytest

DEFAULT_TECHNOLOGIES = set([
    "battery",
    "hydrogen",
    "open_field_pv",
    "wind_onshore_competing",
    "wind_onshore_monopoly",
    "roof_mounted_pv",
    "wind_offshore",
    "hydro_run_of_river",
    "hydro_reservoir",
    "pumped_hydro",
    "demand_elec",
    "nuclear",
])
DIRECTIONAL_PV = set([
    "roof_mounted_pv_s_flat",
    "roof_mounted_pv_n",
    "roof_mounted_pv_e_w",
])
HEAT_TECHS = set([
    "biofuel_boiler",
    "biofuel_tech_heat_to_demand",
    "heat_pump",
    "heat_pump_tech_heat_to_demand",
    "electric_heater",
    "electric_heater_tech_heat_to_demand",
    "hp_heat_storage_small",
    "electric_heater_heat_storage_small",
    "biofuel_heat_storage_small",
    "methane_heat_storage_small",
    "methane_boiler",
    "methane_tech_heat_to_demand",
])
BIOFUEL_TECHS = set([
    "biofuel_supply",
    "electricity_from_biofuel",
])

SYNFUEL_TECHS = set([
    "hydrogen_to_liquids",
    "hydrogen_to_methanol",
    "hydrogen_to_methane",
    "dac",
    "biofuel_to_liquids",
    "biofuel_to_diesel",
    "biofuel_to_methanol",
    "biofuel_to_methane",
    "electrolysis",
    "syn_diesel_converter",
    "syn_methane_converter",
    "syn_kerosene_converter",
    "syn_methanol_converter",
    "syn_diesel_distribution_export",
    "syn_methane_distribution_export",
    "syn_kerosene_distribution_export",
    "syn_methanol_distribution_export",
    "syn_diesel_distribution_import",
    "syn_methane_distribution_import",
    "syn_kerosene_distribution_import",
    "syn_methanol_distribution_import",
])
# Only includes scenarios with non-default technology sets
TECHNOLOGIES = {
    "connected_all_neighbours": DEFAULT_TECHNOLOGIES | set(["ac_transmission"]),
    "connected_entsoe_tyndp": DEFAULT_TECHNOLOGIES | set(["ac_transmission"]),
    "directional-pv": (DEFAULT_TECHNOLOGIES | DIRECTIONAL_PV)
    - set(["roof_mounted_pv"]),
    "shed-load": DEFAULT_TECHNOLOGIES | set(["load_shedding"]),
    "all-overrides": (
        (
            DEFAULT_TECHNOLOGIES
            | DIRECTIONAL_PV
            | set(["load_shedding"])
            | HEAT_TECHS
            | BIOFUEL_TECHS
            | SYNFUEL_TECHS
        )
        - set(["roof_mounted_pv"])
    ),
    "electrified-heat": DEFAULT_TECHNOLOGIES | set(["historic_electrified_heat"]),
    "electrified-biofuel": DEFAULT_TECHNOLOGIES | set(["electrified_biofuel"]),
    "heat": DEFAULT_TECHNOLOGIES | HEAT_TECHS,
    "biofuel": DEFAULT_TECHNOLOGIES | BIOFUEL_TECHS,
    "synfuel": DEFAULT_TECHNOLOGIES | SYNFUEL_TECHS,
}
OPTIONAL_LOCATIONAL_TECHNOLOGIES = ["nuclear"]


@pytest.fixture(scope="function")
def technologies(scenario):
    return TECHNOLOGIES.get(scenario, DEFAULT_TECHNOLOGIES)


def test_model_runs(optimised_model):
    assert optimised_model.results.termination_condition == "optimal"


def test_example_model_runs(optimised_example_model):
    assert optimised_example_model.results.termination_condition == "optimal"


def test_technologies_are_available(
    optimised_model, energy_cap, location, technologies
):
    for technology in technologies:
        if "transmission" in technology:
            assert pd.notna(
                energy_cap.where(energy_cap.techs.str.find(technology) > -1)
                .sum(min_count=1)
                .item()
            )
        elif technology in OPTIONAL_LOCATIONAL_TECHNOLOGIES:
            # don't check the capacity values at each location,
            # since we can't be certain that the technology exists at any specific location
            assert technology in energy_cap.techs
        else:
            assert (technology in energy_cap.techs) and pd.notna(
                energy_cap.sel(locs=location, techs=technology).item()
            )


def test_heat_carrier_exists(model, scenario):
    if scenario in ["heat", "all-overrides"]:
        assert "heat" in model.inputs.carriers
    else:
        assert "heat" not in model.inputs.carriers
