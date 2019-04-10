from pathlib import Path

import pytest
import calliope

PATH_TO_SIMPLE_NATIONAL_MODEL = Path(__file__).parent / "simple-national-model.yaml"
PATH_TO_SIMPLE_REGIONAL_MODEL = Path(__file__).parent / "simple-regional-model.yaml"
PATH_TO_CONNECTED_NATIONAL_MODEL = Path(__file__).parent / "connected-national-model.yaml"
PATH_TO_CONNECTED_REGIONAL_MODEL = Path(__file__).parent / "connected-regional-model.yaml"


@pytest.mark.parametrize("path_to_model", [
    (PATH_TO_SIMPLE_NATIONAL_MODEL), (PATH_TO_SIMPLE_REGIONAL_MODEL)
])
def test_simple_model_runs(path_to_model):
    model = calliope.Model(path_to_model.as_posix())
    model.run()
    assert model.results.termination_condition == "optimal"


@pytest.mark.parametrize("path_to_model", [
    (PATH_TO_CONNECTED_NATIONAL_MODEL), (PATH_TO_CONNECTED_REGIONAL_MODEL)
])
def test_connected_neighbours_model_runs(path_to_model):
    model = calliope.Model(path_to_model.as_posix())
    model.run()
    assert model.results.termination_condition == "optimal"
