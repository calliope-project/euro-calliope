import os
import sys
from pathlib import Path

import calliope
import pandas as pd
import pytest
import xarray as xr
import yaml


def run_test(snakemake):
    with open(
        os.path.join(snakemake.input.test_dir, "..", "resources", "test.yaml")
    ) as f:
        test_config = yaml.safe_load(f)

    override_dict = test_config["test-model"]["overrides"][
        snakemake.wildcards.resolution
    ]
    scenarios = test_config["test-model"]["scenarios"][snakemake.wildcards.resolution]
    subset_time = test_config["test-model"]["subset_time"][
        snakemake.wildcards.resolution
    ]

    exit_code = pytest.main(
        [
            snakemake.input.test_dir,
            f"--html={snakemake.log[0]}",
            "--self-contained-html",
            "--verbose",
        ],
        plugins=[
            _create_config_plugin(snakemake, override_dict, scenarios, subset_time)
        ],
    )
    if exit_code == 0:
        Path(snakemake.output[0]).touch()
    sys.exit(exit_code)


def _create_config_plugin(snakemake, override_dict, scenarios, subset_time):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin:
        @pytest.fixture(scope="session")
        def config(self):
            return snakemake.params.config

        @pytest.fixture(scope="session")
        def scaling_factors(self, config):
            return config["scaling-factors"]

        @pytest.fixture(scope="session")
        def override_dict(self):
            return {"model.subset_time": subset_time, "overrides": override_dict}

        @pytest.fixture(scope="session", params=list(scenarios.keys()))
        def scenario(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def model(self, scenario, override_dict):
            return calliope.Model(
                snakemake.input.example_model,
                scenario=",".join(scenarios[scenario]),
                override_dict=override_dict,
            )

        @pytest.fixture(scope="session")
        def optimised_model(self, model):
            model.run()
            return model

        @pytest.fixture(scope="session")
        def energy_cap(self, optimised_model, scaling_factors):
            return (
                optimised_model.get_formatted_array("energy_cap")
                / scaling_factors["power"]
            )

        @pytest.fixture(
            scope="module",
            params=_read_locs(snakemake.input.capacity_factor_timeseries[0]),
        )
        def location(self, request):
            return request.param

        @pytest.fixture(scope="session")
        def optimised_example_model(self, override_dict):
            model = calliope.Model(
                snakemake.input.example_model,
                override_dict=override_dict,
                scenario=",".join(scenarios["default"]),
            )
            model.run()
            return model

        @pytest.fixture(
            scope="function", params=snakemake.input.capacity_factor_timeseries
        )
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

        @pytest.fixture(scope="module")
        def wind_onshore_capacity_factor_timeseries(self):
            path = self._select_capacity_factor_time_series("wind-onshore")
            return pd.read_csv(path, index_col=0, parse_dates=True)

        def _select_capacity_factor_time_series(self, technology):
            selected = [
                path
                for path in snakemake.input.capacity_factor_timeseries
                if Path(path).name == f"capacityfactors-{technology}.csv"
            ]
            assert len(selected) == 1
            return selected[0]

        @pytest.fixture(scope="module")
        def cop_timeseries(self):
            return pd.read_csv(snakemake.input.cop, index_col=0, parse_dates=True)

        @pytest.fixture(scope="module")
        def unscaled_space_heat_timeseries(self):
            return xr.open_dataset(snakemake.input.unscaled_space_heat)

    return SnakemakeConfigPlugin()


def _read_locs(path_to_cf_timeseries):
    return pd.read_csv(path_to_cf_timeseries, index_col=0, parse_dates=True).columns


if __name__ == "__main__":
    run_test(snakemake)
