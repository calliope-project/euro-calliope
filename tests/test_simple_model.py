from pathlib import Path

import calliope

PATH_TO_SIMPLE_MODEL = Path(__file__).parent / "simple-model.yaml"


def test_simple_model_runs():
    model = calliope.Model(PATH_TO_SIMPLE_MODEL.as_posix())
    model.run()
    assert model.results.termination_condition == "optimal"
