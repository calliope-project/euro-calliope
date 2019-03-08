from pathlib import Path

import pytest
import calliope

PATH_TO_NATIONAL_MODEL = Path(__file__).parent / "simple-national-model.yaml"
PATH_TO_REGIONAL_MODEL = Path(__file__).parent / "simple-regional-model.yaml"


@pytest.mark.parametrize("path_to_model", [
    (PATH_TO_NATIONAL_MODEL), (PATH_TO_REGIONAL_MODEL)
])
def test_simple_model_runs(path_to_model):
    model = calliope.Model(path_to_model.as_posix())
    model.run()
    assert model.results.termination_condition == "optimal"
