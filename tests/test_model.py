def test_model_runs(model):
    model.run()
    assert model.results.termination_condition == "optimal"
