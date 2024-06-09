"""
Use weather data to simulate hourly heat profiles using When2Heat
(https://github.com/oruhnau/when2heat) data processing pipeline.

Functions attributable to When2Heat are explicitly referenced as such in the function docstring.
"""

import os
from typing import Callable, Literal, Union

import numpy as np
import pandas as pd
import xarray as xr

# ASSUME (from When2Heat):
# 1. Below 15 °C, the water heating demand is not defined and assumed to stay constant
HOT_WATER_LOWER_BOUND_TEMP = 15
# 2. the reference temperature is always 30C for hot water demand calcs
HOT_WATER_REF_TEMP = 30
# All locations are separated by the average wind speed with the threshold 4.4 m/s separating windy and not-windy (normal) locations
AVE_WIND_SPEED_THRESHOLD = 4.4


def get_unscaled_heat_profiles(
    path_to_population: str,
    path_to_wind_speed: str,
    path_to_temperature: str,
    path_to_when2heat_params: str,
    first_year: Union[str, int],
    final_year: Union[str, int],
    out_path: str,
) -> None:
    """Produces time series of heat demand profiles with the correct shape, and consistent within themselves, but without meaningful units.

    The profiles need to be scaled so their magnitude matches annual heat demand data in a subsequent step.

    Args:
        path_to_population (str): Gridded population data, which will act as the weighting to got from grid-level profiles to profiles per model unit.
        path_to_wind_speed (str): Gridded wind speed data in m/s.
        path_to_temperature (str): Gridded air temperature data in degrees C.
        path_to_when2heat_params (str): Parameters to convert weather data into demand profiles from the When2Heat project.
        first_year (Union[str, int]): First year of data to include in the profile (inclusive).
        final_year (Union[str, int]): Final year of data to include in the profile (inclusive).
        out_path (str): Path to which data will be saved.
    """
    population = xr.open_dataarray(path_to_population)

    # Weather data is subset by the geographic area covered by model run (given by available population sites)
    temperature_ds = xr.open_dataset(path_to_temperature).sel(site=population.site)
    wind_ds = xr.open_dataset(path_to_wind_speed).sel(site=population.site)

    # Check units
    assert temperature_ds.attrs["unit"].lower() == "degrees c"
    assert wind_ds.attrs["unit"].lower() == "m/s"

    # Only need site-wide mean wind speed for this analysis
    average_wind_speed = wind_ds["wind10m"].mean("time")

    # Subset temperature to the selected year extended by a couple of days either end,
    # so we don't compute values for years we don't need, but keep a buffer for the shifts
    # happening in get_reference_temperature()
    temperature_ds = temperature_ds.sel(
        time=slice(
            str(int(first_year) - 1) + "-12-25", str(int(final_year) + 1) + "-01-05"
        )
    )

    # This is a weighted average temperature from 3 days prior to each day in the timeseries
    # to represent the relative impact of historical daily temperature on the heat demand of each day.
    # See [@Ruhnau:2019] for more information on the method
    reference_temperature = get_reference_temperature(
        temperature_ds["temperature"], time_dim="time"
    )

    # After running get_reference_temperature(), we now subset to get only the target year
    reference_temperature = reference_temperature.sel(
        time=slice(str(first_year), str(final_year))
    )

    # Parameters and how to apply them is based on [@BDEW:2015]
    daily_params = read_daily_parameters(path_to_when2heat_params)
    hourly_params = read_hourly_parameters(path_to_when2heat_params)

    # Get daily demand
    daily_heat = daily(
        reference_temperature, average_wind_speed, daily_params, _heat_function
    )
    daily_hot_water = daily(
        reference_temperature, average_wind_speed, daily_params, _water_function
    )

    # Map profiles to daily demand
    hourly_heat = get_hourly_heat_profiles(
        reference_temperature, daily_heat, hourly_params
    )

    hourly_hot_water = get_hourly_heat_profiles(
        reference_temperature.clip(min=HOT_WATER_REF_TEMP),
        daily_hot_water,
        hourly_params,
    )

    # Space heating demand = total heating demand - hot water demand
    hourly_space = (hourly_heat - hourly_hot_water).clip(min=0)

    # Sanity check that there is more space heating demand in winter than summer
    hourly_space_monthly = hourly_space.groupby("time.month").sum()
    hourly_space_winter = hourly_space_monthly.sel(month=[12, 1, 2]).sum("month")
    hourly_space_summer = hourly_space_monthly.sel(month=[6, 7, 8]).sum("month")
    assert (hourly_space_winter > hourly_space_summer).all()

    # population weighted profiles.
    # NOTE: profile magnitude is now only consistent within each region, not between them
    weight = population / population.sum(["site"])
    grouped_hourly_space = _group_gridcells(hourly_space, weight).rename("space_heat")
    grouped_hourly_water = _group_gridcells(hourly_hot_water, weight).rename(
        "hot_water"
    )

    grouped_hourly_heat = xr.merge([grouped_hourly_space, grouped_hourly_water])
    grouped_hourly_heat.to_netcdf(out_path)


def get_hourly_heat_profiles(
    reference_temperature: xr.DataArray,
    daily_heat: xr.DataArray,
    hourly_params: pd.Series,
) -> xr.DataArray:
    """Convert daily heat demand to hourly profiles.

    Heavily modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py
    to work with xarray datasets and to improve efficiency

    Args:
        reference_temperature (xr.DataArray): Daily reference temperature in degrees C.
        daily_heat (xr.DataArray): Relative daily heat demand per site.
        hourly_params (pd.Series): Parameters from When2Heat to convert from daily to hourly profiles.

    Returns:
        xr.DataArray: Hourly heat demand profiles (which are internally consistent but whose magnitudes are meaningless).
    """

    # get temperature in 5C increments between -15C and +30C
    temperature_increments = (
        (np.ceil((reference_temperature / 5).astype("float64")) * 5)
        .clip(min=-15, max=30)
        .to_series()
    )
    # Profiles are linked to the day's temperature increment, so we align the two
    increment_index = (
        temperature_increments.to_frame("temperature")
        .assign(weekday=temperature_increments.index.get_level_values("time").dayofweek)
        .set_index(["temperature", "weekday"], append=True)
    )
    aligned_increment_to_params = pd.merge(
        hourly_params.rename("param"),
        increment_index,
        left_index=True,
        right_index=True,
    ).squeeze()

    hourly_params_at_all_locations = xr.DataArray.from_series(
        aligned_increment_to_params.droplevel(["weekday", "temperature"])
    )
    # daily heat is multiplied by the hourly parameter value to get the relative heat demand for that hour
    hourly_heat = hourly_params_at_all_locations * daily_heat

    hourly_heat = _hour_and_day_to_datetime(hourly_heat)

    return hourly_heat


def read_daily_parameters(input_path: str) -> pd.DataFrame:
    """Load When2Heat daily parameters.

    Direct copy from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/read.py
    """
    file = os.path.join(input_path, "daily_demand.csv")
    return pd.read_csv(file, sep=";", decimal=",", header=[0, 1], index_col=0)


def read_hourly_parameters(input_path: str) -> pd.DataFrame:
    """Load When2Heat hourly parameters

    Modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/read.py
    to set columns as integer.
    """
    parameters = {}

    parameters["COM"] = _csv_reader("COM", input_path)

    for building_type in ["SFH", "MFH"]:
        parameters[building_type] = (
            _csv_reader(building_type, input_path)
            .rename_axis(index="time")
            .align(parameters["COM"])[0]
        )

    # We postprocess the dataframe from the CSV to make it easier to use in subsequent operations.
    combined_df = pd.concat(
        parameters.values(), keys=parameters.keys(), names=["building"]
    )
    combined_df.columns = combined_df.columns.astype(int).rename("temperature")
    combined_df = combined_df.rename(lambda x: int(x.replace(":00", "")), level="time")
    combined_df = combined_df.rename_axis(index={"time": "hour"})
    return combined_df.stack()


def get_reference_temperature(
    temperature: xr.DataArray, time_dim: str = "time"
) -> xr.DataArray:
    """Get daily reference temperature values which account for the temperature in preceding days using a weighted average.

    Modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py
    to expect xarray not pandas

    Args:
        temperature (xr.DataArray): Hourly temperature in degrees C
        time_dim (str, optional): The name of the hourly time dimension. Defaults to "time".

    Returns:
        xr.DataArray: Daily reference temperatures per site.
    """

    # Daily average
    # pandas manages time resampling much quicker than xarray, so we switch to a dataframe here.
    daily_average = (
        temperature.rename(time=time_dim)
        .to_series()
        .unstack("time")
        .T.resample("1D")
        .mean()
        .stack()
        .to_xarray()
    )

    # Weighted mean, method for which is given in [@Ruhnau:2019]
    return sum(
        (0.5**i) * daily_average.shift({"time": i}).bfill("time") for i in range(4)
    ) / sum(0.5**i for i in range(4))


def daily(
    temperature: xr.DataArray,
    wind: xr.DataArray,
    all_parameters: pd.DataFrame,
    func: Callable,
) -> xr.DataArray:
    """
    Modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py to not prescribe index level names.
    All locations are separated by the average wind speed with the threshold 4.4 m/s
    separating windy and not-windy (normal) locations, then relevant paramaters for that
    windiness are applied in a function derived from [@BDEW:2015] to get "daily heat demand".
    """

    buildings = ["SFH", "MFH", "COM"]

    return xr.concat(
        [
            func(
                temperature.where(wind <= AVE_WIND_SPEED_THRESHOLD),
                all_parameters[(building, "normal")],
            ).fillna(
                func(
                    temperature.where(wind > AVE_WIND_SPEED_THRESHOLD),
                    all_parameters[(building, "windy")],
                )
            )
            for building in buildings
        ],
        dim=pd.Index(buildings, name="building"),
    )


def _group_gridcells(gridded_data: xr.DataArray, weight: xr.DataArray) -> xr.DataArray:
    # `hourly_heat` has dims [x, y, datetime], `weight` has dims [x, y, id],
    # we want a final array with dims [id, datetime]

    return xr.concat(
        [(gridded_data * weight.sel({"id": id})).sum(["site"]) for id in weight.id],
        dim="id",
    )


def _hour_and_day_to_datetime(da: xr.DataArray) -> xr.DataArray:
    "Combine hour and date (a.k.a. 'time') as two dimensions into one datetime ('time') dimension"
    da = da.stack(new_time=["time", "hour"])
    new_time = da.new_time.to_index()
    da.coords["new_time"] = new_time.get_level_values(0) + pd.to_timedelta(
        new_time.get_level_values(1), unit="H"
    )
    return da.rename({"new_time": "time"})


def _csv_reader(
    building_type: Literal["SFH", "MFH", "COM"], input_path: str
) -> pd.DataFrame:
    filename = f"hourly_factors_{building_type}.csv"
    filepath = os.path.join(input_path, filename)

    # MultiIndex for commercial heat because of weekday dependency
    index_col = [0, 1] if building_type == "COM" else 0
    return pd.read_csv(filepath, sep=";", decimal=",", index_col=index_col).apply(
        pd.to_numeric, downcast="float"
    )


def _heat_function(temperature: xr.DataArray, parameters: pd.DataFrame) -> xr.DataArray:
    """A function for the total (space + water) daily heating demand, derived from [@BDEW:2015].

    Direct copy from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py

    Args:
        temperature (xr.DataArray): Reference daily temperature in degrees C.
        parameters (pd.DataFrame): Relevant parameters from When2Heat.

    Returns:
        xr.DataArray: Daily relative heat demand.
    """

    sigmoid = (
        parameters["A"]
        / (1 + (parameters["B"] / (temperature - 40)) ** parameters["C"])
        + parameters["D"]
    )

    linear = xr.concat(
        [parameters[f"m_{i}"] * temperature + parameters[f"b_{i}"] for i in ["s", "w"]],
        dim="param",
    ).max("param")

    return sigmoid + linear


def _water_function(
    temperature: xr.DataArray, parameters: pd.DataFrame
) -> xr.DataArray:
    """A function for the daily water heating demand, derived from [@BDEW:2015].

    Direct copy from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py

    Args:
        temperature (xr.DataArray): Reference daily temperature in degrees C.
        parameters (pd.DataFrame): Relevant parameters from When2Heat.

    Returns:
        xr.DataArray: Daily relative hot water demand.

    """
    celsius_clipped = temperature.clip(min=HOT_WATER_LOWER_BOUND_TEMP)

    return parameters["m_w"] * celsius_clipped + parameters["b_w"] + parameters["D"]


if __name__ == "__main__":
    get_unscaled_heat_profiles(
        path_to_population=snakemake.input.population,
        path_to_wind_speed=snakemake.input.wind_speed,
        path_to_temperature=snakemake.input.temperature,
        path_to_when2heat_params=snakemake.input.when2heat,
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        out_path=snakemake.output[0],
    )
