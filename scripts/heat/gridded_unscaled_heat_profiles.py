"""
Use weather data to simulate hourly heat profiles using When2Heat
(https://github.com/oruhnau/when2heat) data processing pipeline.
Functions attributable to When2Heat are explictily referenced as such
in the function docstring.
"""

import os
from typing import Callable, Union

import numpy as np
import pandas as pd
import xarray as xr


def get_unscaled_heat_profiles(
    path_to_population: str,
    path_to_wind_speed: str,
    path_to_temperature: str,
    path_to_when2heat_params: str,
    lat_name: str,
    lon_name: str,
    out_path: str,
    year: Union[str, int],
) -> None:
    population = xr.open_dataarray(path_to_population)
    temperature_ds = xr.open_dataset(path_to_temperature).rename({
        lon_name: "x",
        lat_name: "y",
    })
    # Only need site-wide mean wind speed for this analysis
    wind_ds = xr.open_dataset(path_to_wind_speed).rename({lon_name: "x", lat_name: "y"})
    average_wind_speed = wind_ds["wind_speed"].mean("time")

    # Subset temperature to the selected year extended by a couple of days either end,
    # so we don't compute values for years we don't need, but keep a buffer for the shifts
    # happening in get_reference_temperature()
    temperature_ds = temperature_ds.sel(
        time=slice(str(int(year) - 1) + "-12-25", str(int(year) + 1) + "-01-05")
    )

    # This is a weighted average temperature from 3 days prior to each day in the timeseries
    # to represent the relative impact of historical daily temperature on the heat demand of each day.
    # See [@Ruhnau:2019] for more information on the method
    reference_temperature = get_reference_temperature(
        temperature_ds["temperature"], time_dim="time"
    )

    # After running get_reference_temperature(), we now subset to get only the target year
    reference_temperature = reference_temperature.sel(time=str(year))

    # Parameters and how to apply them is based on [@BDEW:2015]
    daily_params = read_daily_parameters(path_to_when2heat_params)
    hourly_params = read_hourly_parameters(path_to_when2heat_params)

    # Get daily demand
    daily_heat = get_daily_heat_demand(
        reference_temperature, average_wind_speed, daily_params
    )

    # Map profiles to daily demand
    hourly_heat = get_hourly_heat_profiles(
        reference_temperature, daily_heat, hourly_params
    )

    # population weighted profiles.
    # NOTE: profile magnitude is now only consistent within each region, not between them
    weight = population / population.sum(["site"])
    # `hourly_heat` has dims [x, y, datetime], `weight` has dims [x, y, id],
    # we want a final array with dims [id, datetime]
    grouped_hourly_heat = xr.concat(
        [(hourly_heat * weight.sel({"id": id})).sum(["site"]) for id in weight.id],
        dim="id",
    )

    grouped_hourly_heat.to_netcdf(out_path)


def get_hourly_heat_profiles(
    reference_temperature: Union[int, float],
    daily_heat: xr.DataArray,
    hourly_params: pd.DataFrame,
) -> xr.DataArray:
    """
    reference_temperature: temperature in degrees C

    """
    # Heavily modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py
    # to work with xarray datasets and to improve efficiency

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
    aligned_increment_to_params = hourly_params.align(increment_index)[0]
    hourly_params_at_all_locations = xr.DataArray.from_series(
        aligned_increment_to_params.droplevel(["weekday", "temperature"])
    )
    # daily heat is multiplied by the hourly paramater value to get the relative heat demand for that hour
    hourly_heat = hourly_params_at_all_locations * daily_heat

    hourly_heat = hour_and_day_to_datetime(hourly_heat)

    return hourly_heat


def hour_and_day_to_datetime(da: xr.DataArray) -> xr.DataArray:
    da = da.stack(new_time=["time", "hour"])
    new_time = da.new_time.to_index()
    da.coords["new_time"] = new_time.get_level_values(0) + pd.to_timedelta(
        new_time.get_level_values(1), unit="H"
    )
    return da.rename({"new_time": "time"})


def read_daily_parameters(input_path: str) -> pd.DataFrame:
    # Direct copy from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/read.py
    file = os.path.join(input_path, "daily_demand.csv")
    return pd.read_csv(file, sep=";", decimal=",", header=[0, 1], index_col=0)


def read_hourly_parameters(input_path: str) -> pd.DataFrame:
    # Modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/read.py
    # to set columns as integer
    parameters = {}

    def _csv_reader(building_type):
        filename = f"hourly_factors_{building_type}.csv"
        filepath = os.path.join(input_path, filename)

        # MultiIndex for commercial heat because of weekday dependency
        index_col = [0, 1] if building_type == "COM" else 0
        return pd.read_csv(filepath, sep=";", decimal=",", index_col=index_col).apply(
            pd.to_numeric, downcast="float"
        )

    parameters["COM"] = _csv_reader("COM")

    for building_type in ["SFH", "MFH"]:
        parameters[building_type] = (
            _csv_reader(building_type)
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


def get_reference_temperature(temperature: xr.DataArray, time_dim: str = "time"):
    # Modified from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py
    # to expect xarray not pandas

    # Daily average
    daily_average = temperature.resample(time="1D").mean(time_dim)

    # Weighted mean, method for which is given in [@Ruhnau:2019]
    return sum(
        (0.5**i) * daily_average.shift({time_dim: i}).bfill(time_dim) for i in range(4)
    ) / sum(0.5**i for i in range(4))


def get_daily_heat_demand(
    temperature: xr.DataArray, average_wind: xr.DataArray, all_parameters: pd.DataFrame
) -> xr.DataArray:
    # Direct copy from https://github.com/oruhnau/when2heat/blob/351bd1a2f9392ed50a7bdb732a103c9327c51846/scripts/demand.py

    def heat_function(t, parameters):
        # BDEW et al. 2015 describes this function combining parameters
        celsius = t - 273.15  # The temperature input is in Kelvin

        sigmoid = (
            parameters["A"]
            / (1 + (parameters["B"] / (celsius - 40)) ** parameters["C"])
            + parameters["D"]
        )

        linear = xr.concat(
            [parameters[f"m_{i}"] * celsius + parameters[f"b_{i}"] for i in ["s", "w"]],
            dim="param",
        ).max("param")

        return sigmoid + linear

    return daily(temperature, average_wind, all_parameters, heat_function)


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
                temperature.where(wind <= 4.4), all_parameters[(building, "normal")]
            ).fillna(
                func(temperature.where(wind > 4.4), all_parameters[(building, "windy")])
            )
            for building in buildings
        ],
        dim=pd.Index(buildings, name="building"),
    )


if __name__ == "__main__":
    get_unscaled_heat_profiles(
        path_to_population=snakemake.input.population,
        path_to_wind_speed=snakemake.input.wind_speed,
        path_to_temperature=snakemake.input.temperature,
        path_to_when2heat_params=snakemake.input.when2heat,
        lat_name=snakemake.params.lat_name,
        lon_name=snakemake.params.lon_name,
        year=snakemake.wildcards.year,
        out_path=snakemake.output[0],
    )
