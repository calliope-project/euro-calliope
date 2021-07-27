from pathlib import Path
import sys
import yaml

import pytest
import calliope
import pandas as pd


def run_test(path_to_test_dir, path_to_output, path_to_model, path_to_example_model, paths_to_cf_timeseries,
             path_to_annual_energy_balances, config):
    exit_code = pytest.main(
        [
            path_to_test_dir,
            f"--html={path_to_output}",
            "--self-contained-html",
        ],
        plugins=[
            _create_config_plugin(
                path_to_model=path_to_model,
                path_to_example_model=path_to_example_model,
                paths_to_cf_timeseries=paths_to_cf_timeseries,
                path_to_annual_energy_balances=path_to_annual_energy_balances,
                config=config
            )
        ]
    )
    sys.exit(exit_code)


def _create_config_plugin(path_to_model, path_to_example_model, paths_to_cf_timeseries, path_to_annual_energy_balances, config):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="session")
        def config(self):
            return dict(config)

        @pytest.fixture(scope="session")
        def scaling_factors(self, config):
            return config["scaling-factors"]

        @pytest.fixture(scope="session", params=_read_scenario_names_from_yaml(path_to_model))
        def scenario(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def model(self, scenario):
            return calliope.Model(path_to_model, scenario=scenario)

        @pytest.fixture(scope="session")
        def optimised_model(self, model):
            model.run()
            return model

        @pytest.fixture(scope="session")
        def energy_cap(self, optimised_model, scaling_factors):
            return optimised_model.get_formatted_array("energy_cap") / scaling_factors["power"]

        @pytest.fixture(scope="module", params=_read_locs(paths_to_cf_timeseries[0]))
        def location(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def optimised_example_model(self):
            model = calliope.Model(path_to_example_model)
            model.run()
            return model

        @pytest.fixture(scope="function", params=paths_to_cf_timeseries)
        def capacity_factor_timeseries(self, request):
            return pd.read_csv(request.param, index_col=0, parse_dates=True)

        @pytest.fixture(scope="module")
        def open_field_pv_capacity_factor_timeseries(self):
            path = self._select_capacity_factor_time_series("open-field-pv")
            return pd.read_csv(path, index_col=0, parse_dates=True)

        @pytest.fixture(scope="module")
        def rooftop_pv_capacity_factor_timeseries(self):
            path = self._select_capacity_factor_time_series("rooftop-pv")
            return pd.read_csv(path, index_col=0, parse_dates=True)

        def _select_capacity_factor_time_series(self, technology):
            selected = [path for path in paths_to_cf_timeseries
                        if Path(path).name == f"capacityfactors-{technology}.csv"]
            assert len(selected) == 1
            return selected[0]

        @pytest.fixture(scope="module", params=path_to_annual_energy_balances)
        def annual_energy_balances(self, request):
            return pd.read_csv(request.param, index_col=0, parse_dates=True)

    return SnakemakeConfigPlugin()


def _read_scenario_names_from_yaml(path_to_model):
    with open(path_to_model, 'r') as stream:
        model = yaml.safe_load(stream)
    return model["scenarios"].keys()


def _read_locs(path_to_cf_timeseries):
    return pd.read_csv(path_to_cf_timeseries, index_col=0, parse_dates=True).columns


if __name__ == "__main__":
    run_test(
        path_to_test_dir=snakemake.input.test_dir,
        path_to_model=snakemake.input.model,
        path_to_example_model=snakemake.input.example_model,
        paths_to_cf_timeseries=snakemake.input.capacity_factor_timeseries,
        path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
        config=snakemake.params.config,
        path_to_output=snakemake.output[0]
    )
