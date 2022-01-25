import pandas as pd
import pytest

DEFAULT_TECHNOLOGIES = set([
    "battery", "hydrogen", "open_field_pv", "wind_onshore_competing", "wind_onshore_monopoly",
    "roof_mounted_pv", "wind_offshore", "hydro_run_of_river", "hydro_reservoir", "pumped_hydro",
    "biofuel", "demand_elec"
])
DIRECTIONAL_PV = set(["roof_mounted_pv_s_flat", "roof_mounted_pv_n", "roof_mounted_pv_e_w"])

# Only includes scenarios with non-default technology sets
TECHNOLOGIES = {
    "connected_all_neighbours": DEFAULT_TECHNOLOGIES | set(["ac_transmission"]),
    "connected_entsoe_tyndp": DEFAULT_TECHNOLOGIES | set(["ac_transmission"]),
    "directional-pv": (DEFAULT_TECHNOLOGIES | DIRECTIONAL_PV) - set(["roof_mounted_pv"]),
    "shed-load": DEFAULT_TECHNOLOGIES | set(["load_shedding"]),
    "all-overrides": (
        (DEFAULT_TECHNOLOGIES | DIRECTIONAL_PV | set(["load_shedding"])) - set(["roof_mounted_pv"])
    ),
}


@pytest.fixture(scope="function")
def technologies(scenario):
    return TECHNOLOGIES.get(scenario, DEFAULT_TECHNOLOGIES)


def test_model_runs(optimised_model):
    assert optimised_model.results.termination_condition == "optimal"


def test_example_model_runs(optimised_example_model):
    assert optimised_example_model.results.termination_condition == "optimal"


def test_technologies_are_available(energy_cap, location, technologies):
    for technology in technologies:
        if "transmission" in technology:
            assert pd.notna(
                energy_cap.where(energy_cap.techs.str.find(technology) > -1).sum(min_count=1).item()
            )
        else:
            assert (
                (technology in energy_cap.techs)
                and pd.notna(energy_cap.sel(locs=location, techs=technology).item())
            )
