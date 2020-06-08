import sys

import pytest
import calliope
import pandas as pd


def run_test(path_to_output, path_to_connected_model, path_to_disconnected_model,
             paths_to_cf_timeseries, config):
    exit_code = pytest.main(
        [
            f"--html={path_to_output}",
            f"--self-contained-html",
        ],
        plugins=[
            _create_config_plugin(
                path_to_connected_model=path_to_connected_model,
                path_to_disconnected_model=path_to_disconnected_model,
                paths_to_cf_timeseries=paths_to_cf_timeseries,
                config=config
            )
        ]
    )
    sys.exit(exit_code)


def _create_config_plugin(path_to_connected_model, path_to_disconnected_model,
                          paths_to_cf_timeseries, config):
    """Creates fixtures from Snakemake configuration."""

    class SnakemakeConfigPlugin():

        @pytest.fixture(scope="function")
        def config(self):
            return dict(config)

        @pytest.fixture(scope="function", params=[path_to_connected_model, path_to_disconnected_model])
        def model(self, request):
            return calliope.Model(request.param)

        @pytest.fixture(scope="function", params=paths_to_cf_timeseries)
        def capacity_factor_timeseries(self, request):
            return pd.read_csv(request.param, index_col=0, parse_dates=True)

    return SnakemakeConfigPlugin()


if __name__ == "__main__":
    run_test(
        path_to_connected_model=snakemake.input.connected_model,
        path_to_disconnected_model=snakemake.input.disconnected_model,
        paths_to_cf_timeseries=snakemake.input.capacity_factor_timeseries,
        config=snakemake.params.config,
        path_to_output=snakemake.output[0]
    )
