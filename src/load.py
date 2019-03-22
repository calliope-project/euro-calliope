"""Generate Calliope load time series."""
import pandas as pd
import geopandas as gpd


def electricity_demand(path_to_units, path_to_electricity_load, scaling_factor, path_to_result):
    """Generate load time series for every location."""
    units = gpd.read_file(path_to_units)
    national_load = pd.read_csv(path_to_electricity_load, index_col=0)

    units["fraction_of_national_load"] = units.groupby("country_code").demand_twh_per_year.transform(
        lambda x: x / x.sum()
    )

    pd.concat(
        [unit_time_series(unit, national_load, scaling_factor) for _, unit in units.iterrows()],
        axis=1
    ).to_csv(path_to_result)


def unit_time_series(unit, national_load, scaling_factor):
    country_code = unit.country_code
    multiplier = unit.fraction_of_national_load
    unit_ts = national_load.loc[:, country_code].copy() * multiplier * (-1) * scaling_factor
    unit_ts.name = unit.id.replace(".", "-")
    return unit_ts


if __name__ == '__main__':
    electricity_demand(
        path_to_units=snakemake.input.units,
        path_to_electricity_load=snakemake.input.national_load,
        path_to_result=snakemake.output[0],
        scaling_factor=snakemake.params.scaling_factor
    )
