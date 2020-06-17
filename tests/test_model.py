def test_model_runs(model):
    model.run()
    assert model.results.termination_condition == "optimal"


def test_example_model_runs(example_model):
    example_model.run()
    assert example_model.results.termination_condition == "optimal"
