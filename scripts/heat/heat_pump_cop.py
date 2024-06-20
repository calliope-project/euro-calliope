from typing import Union

import pandas as pd
import xarray as xr
from eurocalliopelib import utils

HEAT_PUMP_TYPE_SOURCE = {"ashp": "air", "gshp": "ground"}


def cop(
    path_to_temperature_air: str,
    path_to_temperature_ground: str,
    path_to_population: str,
    path_to_heat_pump_characteristics: str,
    sink_temperature: dict[str, int],
    space_heat_sink_shares: dict[str, float],
    correction_factor: float,
    heat_pump_shares: dict[str, float],
    first_year: Union[int, str],
    final_year: Union[int, str],
    path_to_output: str,
):
    """Calculate heat pump Coefficient of Performance (COP) based on manufacturer data.

    COP is calculated for air-source and ground-source heat pumps according to different source temperature data.

    COP is aggregated from grid-cells to model resolution using a population-weighted sum.

    Args:
        path_to_temperature_air (str): Gridded air temperature timeseries data.
        path_to_temperature_ground (str): Gridded ground/soil temperature timeseries data.
        path_to_population (str):  Gridded population data with `id` dimension that defined model resolution unit IDs.
        path_to_heat_pump_characteristics (str):  Manufacturer data on heat pump characteristics across a product range.
        sink_temperature (dict[str, int]): Working temperature for different heating methods (the temperature 'sink' of a heat pump).
        space_heat_sink_shares (dict[str, float]): Share of different space heating methods assumed for the building stock.
        correction_factor (float): Factor with which to downrate heat pump performance to go from manufacturer data to "operational" performance.
        heat_pump_shares (dict[str, float]): Share of air- vs ground-source heat pumps in the market.
        first_year (Union[int, str]): First year of data to include in the profile (inclusive).
        final_year (Union[int, str]): Final year of data to include in the profile (inclusive).
        path_to_output (str): Output to which COP timeseries data will be saved.
    """
    # Initial fast-fail checks.
    assert (
        sum(heat_pump_shares.values()) == 1
    ), "Heat pump technology shares must add up to 1."
    assert (
        sum(space_heat_sink_shares.values()) == 1
    ), "Space heating sink method shares must add up to 1."

    population = xr.open_dataarray(path_to_population)

    temperature_ds = xr.merge(
        _load_temperature_data(filepath, first_year, final_year)
        for filepath in [path_to_temperature_air, path_to_temperature_ground]
    )

    # 1. Get characteristics per sink method
    pre_grouped_heat_pump_characteristics = (
        xr.open_dataarray(path_to_heat_pump_characteristics)
        .mean("product")  # ASSUME: take the average of all heat pump products
        .sel(data_type="COP")
        .interp(sink_temp=list(sink_temperature.values()))
        .assign_coords(sink_temp=list(sink_temperature.keys()))
    )
    # 2. Combine sink methods into space heating and hot water end uses,
    # using weightings for space heating (hot water is a distinct sink method already)
    sink_method_shares = (
        pd.Series({"hot-water": 1, **space_heat_sink_shares})
        .rename_axis(index="sink_temp")
        .to_xarray()
    )
    space_heat_renamer = {k: "space_heat" for k in space_heat_sink_shares}
    heat_pump_characteristics = utils.rename_and_groupby(
        pre_grouped_heat_pump_characteristics * sink_method_shares,
        {"hot-water": "hot_water", **space_heat_renamer},
        "sink_temp",
        "end_use",
        dropna=False,
    )

    cop_ashp = temperature_to_cop(
        heat_pump_characteristics.sel(source="air"),
        temperature_ds["temperature"],
        correction_factor,
    )

    cop_gshp = temperature_to_cop(
        heat_pump_characteristics.sel(source="ground"),
        # ASSUME: 5C decrease to account for soil to brine heat transfer
        temperature_ds["tsoil5"] - 5,
        correction_factor,
    )

    cop = cop_ashp * heat_pump_shares["ashp"] + cop_gshp * heat_pump_shares["gshp"]
    # We infill with ASHP COP for gridcells that have no GSHP data.
    # These tend to be gridcells covering areas with no/limited land.
    cop = cop.fillna(cop_ashp)

    # Sanity check that there is a. higher COP in summer than winter, b. no COP < 1 (worse than direct electrical heating)
    cop_monthly = cop.groupby("time.month").mean()
    cop_winter = cop_monthly.sel(month=[12, 1, 2]).mean("month")
    cop_summer = cop_monthly.sel(month=[6, 7, 8]).mean("month")
    assert (
        cop_summer > cop_winter
    ).all(), "Found higher heating COP values in winter than in summer."
    assert (cop >= 1).all(), "Found improbably low heat pump COP values (< 1)."

    # population weighted COP.
    weight = population / population.sum(["site"])
    # `ds_cop` has dims [site, time], `weight` has dims [site, id], we want a final array with dims [id, time]
    ds_cop_grouped = xr.concat(
        [(cop * weight.sel({"id": id})).sum(["site"]) for id in weight.id],
        dim="id",
    )

    ds_cop_grouped.to_dataset(dim="end_use").to_netcdf(path_to_output)


def temperature_to_cop(
    heat_pump_characteristics: xr.DataArray,
    temperature_celsius: xr.DataArray,
    correction_factor: float,
) -> xr.DataArray:
    """
    Interpolate heat pump temperature-COP relationship to the gridded weather temperature profiles.
    """

    # The range of source temperatures covered in the characteristic data depends on the heat pump type
    source_cop = correction_factor * heat_pump_characteristics.dropna("source_temp")

    return source_cop.interp(
        {"source_temp": temperature_celsius}, kwargs={"fill_value": "extrapolate"}
    )


def _load_temperature_data(
    path_to_temperature_data: str,
    first_year: Union[int, str],
    final_year: Union[int, str],
) -> xr.Dataset:
    "Load xarray dataset and check that units are in the correct unit"
    ds = xr.open_dataset(path_to_temperature_data).sel(
        time=slice(str(first_year), str(final_year))
    )
    assert ds.attrs["unit"].lower() == "degrees c"
    return ds


if __name__ == "__main__":
    cop(
        path_to_temperature_air=snakemake.input.temperature_air,
        path_to_temperature_ground=snakemake.input.temperature_ground,
        path_to_population=snakemake.input.population,
        path_to_heat_pump_characteristics=snakemake.input.heat_pump_characteristics,
        sink_temperature=snakemake.params.sink_temperature,
        space_heat_sink_shares=snakemake.params.space_heat_sink_shares,
        correction_factor=snakemake.params.correction_factor,
        heat_pump_shares=snakemake.params.heat_pump_shares,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        path_to_output=snakemake.output[0],
    )
