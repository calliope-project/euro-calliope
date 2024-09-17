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
    "heat_pump",
    "electric_heater",
    "heat_pump_storage",
    "electric_heater_storage",
    "biofuel_boiler_storage",
    "methane_boiler_storage",
    "methane_boiler",
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
])
# Only includes scenarios with non-default technology sets
TECHNOLOGIES = {
    # TODO: work out way to check for transmission techs
    "connected_all_neighbours": DEFAULT_TECHNOLOGIES,
    "connected_entsoe_tyndp": DEFAULT_TECHNOLOGIES,
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


def test_technologies_are_available(flow_cap, location, technologies):
    for technology in technologies:
        if "link_" in technology:
            assert pd.notna(flow_cap.sel(techs=technology).sum(min_count=1).item())
        elif technology in OPTIONAL_LOCATIONAL_TECHNOLOGIES:
            # don't check the capacity values at each location,
            # since we can't be certain that the technology exists at any specific location
            assert technology in flow_cap.techs
        else:
            assert (technology in flow_cap.techs) and pd.notna(
                flow_cap.sel(nodes=location, techs=technology)
                .sum("carriers", min_count=1)
                .item()
            )


def test_heat_carrier_exists(model, scenario):
    if scenario in ["heat", "all-overrides"]:
        assert "heat" in model.inputs.carriers
    else:
        assert "heat" not in model.inputs.carriers
