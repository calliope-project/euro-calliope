import pandas as pd
import xarray as xr
from eurocalliopelib import utils

# We define independent COP bounds for data checking in case errors creep in with the loaded COP data
MIN_COP = 1
MAX_COP = 6


def create_heat_demand_timeseries(
    path_to_annual_demand: str,
    path_to_unscaled_demand_profiles: str,
    path_to_heat_pump_cop: str,
    power_scaling_factor: float,
    sfh_mfh_ratio: dict,
    historic: bool,
    path_to_output: str,
) -> None:
    """Create a timeseries of electricity demand to meet building space heat and hot water demands.

    ASSUME: all demand is met by heat pumps without any storage buffer.
    ASSUME: if calculating historic electrical heating requirements, annual electricity demand is distributed using the demand profile, ignoring heat pump COP.

    Args:
        path_to_annual_demand (str):
            Annual heat demand by year, resolution ID, end-use, and building category.
        path_to_unscaled_demand_profiles (str):
            Hourly heat demand profiles per year by resolution ID, end-use, and building category.
            Profiles will be normalised before scaling with annual demand, so their absolute magnitudes will be ignored.
        path_to_heat_pump_cop (str):
            Heat pump Coefficient of Performance timeseries per year by resolution ID and end-use.
        power_scaling_factor (float): _description_
            Workflow electricity units scaling factor.
        sfh_mfh_ratio (dict):
            Ratio of single- and multi-family households, used to combine respective unscaled demand profiles into one "household" profile.
        historic (bool):
            If True, calculate electricity demand timeseries according to historic electricity consumption to meet heat demand, not historic heat demand.
        path_to_output (str): Path to which heat demand profiles will be saved.
    """
    assert (
        sum(sfh_mfh_ratio.values()) == 1
    ), "Household type (single- vs multi-family home) ratios must add up to 1."

    annual_demand = (
        pd.read_csv(path_to_annual_demand, index_col=[0, 1, 2])
        .rename_axis(columns="id")
        .stack()
        .unstack("end_use")
        .to_xarray()[["space_heat", "hot_water"]]
    )
    # Read annual heat demand as an xarray dataset
    unscaled_demand_profiles = xr.open_dataset(path_to_unscaled_demand_profiles)

    sfh_mfh_ratio_da = (
        pd.Series({"COM": 1, **sfh_mfh_ratio}).rename_axis(index="building").to_xarray()
    )
    household_building_renamer = {k: "household" for k in sfh_mfh_ratio}
    grouped_unscaled_demand = utils.rename_and_groupby(
        unscaled_demand_profiles * sfh_mfh_ratio_da,
        {"COM": "commercial", **household_building_renamer},
        "building",
        "cat_name",
    )
    scaled_demand_profiles = grouped_unscaled_demand.groupby("time.year").apply(
        _scale_demand, annual_demand=annual_demand
    )

    if historic:
        electricity_demand_twh = scaled_demand_profiles
    else:
        heat_pump_cop = xr.open_dataset(path_to_heat_pump_cop)
        electricity_demand_twh = -1 * scaled_demand_profiles / heat_pump_cop

    # After combining demand and COP, we can aggregate space heating and hot water profiles into one.
    electricity_demand_twh = electricity_demand_twh.to_array("end_use").sum("end_use")
    _check_data(annual_demand, electricity_demand_twh)
    # TWh -> MWh, then scale
    electricity_demand_scaled = electricity_demand_twh * 1e6 * power_scaling_factor

    df_timeseries = electricity_demand_scaled.to_series().unstack("id")

    assert not df_timeseries.isna().any(
        axis=None
    ), "There are NaN values in the timeseries dataframe"

    df_timeseries.to_csv(path_to_output)


def _scale_demand(
    one_year_profile: xr.Dataset, annual_demand: xr.Dataset
) -> xr.Dataset:
    "Scale demand in each year and sum different building types into one profile."
    year = one_year_profile.time.dt.year[0]
    normalised_profile = one_year_profile / one_year_profile.sum("time")
    demand = normalised_profile * annual_demand.sel(year=year).drop("year")
    return demand.sum("cat_name")


def _check_data(heat_demand_twh: xr.Dataset, electricity_demand_twh: xr.DataArray):
    "Check data after completing all processing steps"
    annual_heat_demand = heat_demand_twh.to_array("end_use").sum([
        "cat_name",
        "end_use",
    ])

    cop = annual_heat_demand / abs(electricity_demand_twh).groupby("time.year").sum()

    assert ((cop.round(1) >= MIN_COP) & (cop.round(1) <= MAX_COP)).all()


if __name__ == "__main__":
    create_heat_demand_timeseries(
        path_to_annual_demand=snakemake.input.annual_demand,
        path_to_unscaled_demand_profiles=snakemake.input.unscaled_demand_profiles,
        path_to_heat_pump_cop=snakemake.input.heat_pump_cop,
        power_scaling_factor=snakemake.params.power_scaling_factor,
        sfh_mfh_ratio=snakemake.params.sfh_mfh_ratio,
        historic=snakemake.params.historic,
        path_to_output=snakemake.output[0],
    )
