from pathlib import Path
from typing import Literal, Union

import pandas as pd
import xarray as xr
from eurocalliopelib import utils

HEAT_PUMP_TYPE_SOURCE = {"ashp": "air", "gshp": "ground"}


def cop(
    path_to_temperature: str,
    path_to_population: str,
    path_to_heat_pump_characteristics: str,
    sink_temperature: dict,
    space_heat_sink_ratio: dict,
    correction_factor: float,
    lat_name: str,
    lon_name: str,
    heat_pump_type: Literal["ashp", "gshp"],
    year: Union[str, int],
    path_to_output: str,
):
    """Calculate heat pump Coefficient of Performance (COP) based on manufacturer data.

    COP is calculated for air-source and ground-source heat pumps according to different source temperature data.

    COP is aggregated from grid-cells to model resolution using a population-weighted sum.

    Args:
        path_to_temperature (str): Gridded temperature timeseries data.
        path_to_population (str): Gridded population data with `id` dimension that defined model resolution unit IDs.
        path_to_heat_pump_characteristics (str): Manufacturer data on heat pump characteristics across a product range.
        sink_temperature (dict): Working temperature for different heating methods (the temperature 'sink' of a heat pump).
        space_heat_sink_ratio (dict): Ratio of different space heating methods assumed for the building stock.
        correction_factor: Factor with which to downrate heat pump performance to go from manufacturer data to "operational" performance.
        lat_name (str): Name of the latitude dimension in the gridded datasets.
        lon_name (str):  Name of the longitude dimension in the gridded datasets.
        heat_pump_type (Literal[ashp, gshp]): Heat pump type being modelled.
        year (Union[str, int]): Timeseries year to select.
        path_to_output (str): Output to which COP timeseries data will be saved.
    """
    population = xr.open_dataarray(path_to_population)
    temperature_ds = (
        xr.open_dataset(path_to_temperature)
        .rename({lon_name: "x", lat_name: "y"})
        .sel(time=str(year))
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
    assert (
        sum(space_heat_sink_ratio.values()) == 1
    ), "Space heating sink method ratios must add up to 1."
    sink_method_ratios = (
        pd.Series({"hot-water": 1, **space_heat_sink_ratio})
        .rename_axis(index="sink_temp")
        .to_xarray()
    )
    space_heat_renamer = {k: "space_heat" for k in space_heat_sink_ratio}
    heat_pump_characteristics = utils.rename_and_groupby(
        pre_grouped_heat_pump_characteristics * sink_method_ratios,
        {"hot-water": "hot_water", **space_heat_renamer},
        "sink_temp",
        "end_use",
        dropna=False,
    )

    ds_cop = temperature_to_cop(
        heat_pump_characteristics,
        temperature_ds,
        correction_factor,
        HEAT_PUMP_TYPE_SOURCE[heat_pump_type],
        Path(path_to_temperature).stem,
    )

    # Sanity check that there is higher COP in summer than winter
    cop_monthly = ds_cop.groupby("time.month").sum()
    cop_winter = cop_monthly.sel(month=[12, 1, 2]).mean("month")
    cop_summer = cop_monthly.sel(month=[6, 7, 8]).mean("month")
    assert (cop_summer > cop_winter).all()

    # population weighted COP.
    weight = population / population.sum(["site"])
    # `ds_cop` has dims [site, time], `weight` has dims [site, id],
    # we want a final array with dims [id, time]
    ds_cop_grouped = xr.concat(
        [(ds_cop * weight.sel({"id": id})).sum(["site"]) for id in weight.id],
        dim="id",
    )

    # You can never have a COP < 1 (direct electrical heating)
    assert (
        ds_cop_grouped >= 1
    ).all(), "Improbably low heat pump COP values (< 1) found."

    ds_cop_grouped.to_dataset(dim="end_use").to_netcdf(path_to_output)


def temperature_to_cop(
    heat_pump_characteristics: xr.DataArray,
    temperature_ds: xr.Dataset,
    correction_factor: float,
    source: Literal["ground", "air"],
    weather_var: str,
):
    """
    Interpolate heat pump temperature-COP relationship to the gridded weather temperature profiles.
    """
    if source == "air":
        temperature = temperature_ds[weather_var]
    elif source == "ground":
        # ASSUME: 5C decrease to account for soil to brine heat transfer
        temperature = temperature_ds[weather_var] - 5

    temperature_c = temperature - 273.15 if (temperature > 100).all() else temperature

    # The range of source temperatures covered in the characteristic data depends on the heat pump type
    source_cop = correction_factor * heat_pump_characteristics.dropna(
        "source_temp"
    ).sel(source=source)

    return source_cop.interp(
        {"source_temp": temperature_c}, kwargs={"fill_value": "extrapolate"}
    )


if __name__ == "__main__":
    cop(
        path_to_temperature=snakemake.input.temperature,
        path_to_population=snakemake.input.population,
        path_to_heat_pump_characteristics=snakemake.input.heat_pump_characteristics,
        sink_temperature=snakemake.params.sink_temperature,
        space_heat_sink_ratio=snakemake.params.space_heat_sink_ratio,
        correction_factor=snakemake.params.correction_factor,
        lat_name=snakemake.params.lat_name,
        lon_name=snakemake.params.lon_name,
        heat_pump_type=snakemake.wildcards.heat_pump_type,
        year=snakemake.wildcards.year,
        path_to_output=snakemake.output[0],
    )
