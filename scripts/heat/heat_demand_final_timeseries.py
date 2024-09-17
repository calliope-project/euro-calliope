import numpy as np
import pandas as pd
import xarray as xr
from eurocalliopelib import utils


def scale_heat_demand_profiles(
    annual_demand_twh: xr.Dataset,
    unscaled_demand_profiles: xr.Dataset,
    sfh_mfh_shares: dict,
    model_scaling_factor: float,
) -> xr.DataArray:
    """Create a timeseries of demands for space heat and hot water demands across all building types.

    ASSUME: if calculating historic electrical heating requirements, annual electricity demand is distributed using the demand profile.
    It therefore ignores heat pump COP profiles and the possibility of a storage buffer.

    Args:
        annual_demand_twh (xr.Dataset):
            Annual heat demand in TWh by year, resolution ID, end-use, and building category.
        unscaled_demand_profiles (xr.Dataset):
            Hourly heat demand profiles per year by resolution ID, end-use, and building category.
            Profiles will be normalised before scaling with annual demand, so their absolute magnitudes will be ignored.
        sfh_mfh_shares (dict):
            Share of single- and multi-family households, used to combine respective unscaled demand profiles into one "household" profile.
        model_scaling_factor (float):
            Scaling factor to go from MWh to the units of energy used in the final Calliope model.
    Returns:
        xr.DataArray: `unscaled_demand_profiles` merged across building types and end uses, and scaled to have an annual sum equal to `annual_demand`.
    """
    assert np.isclose(
        sum(sfh_mfh_shares.values()), 1
    ), "Household type (single- vs multi-family home) shares must add up to 1."

    sfh_mfh_shares_da = (
        pd.Series({"COM": 1, **sfh_mfh_shares})
        .rename_axis(index="building")
        .to_xarray()
    )
    household_building_renamer = {k: "household" for k in sfh_mfh_shares}
    grouped_unscaled_demand = utils.rename_and_groupby(
        unscaled_demand_profiles * sfh_mfh_shares_da,
        {"COM": "commercial", **household_building_renamer},
        "building",
        "cat_name",
    )
    scaled_demand_profiles = grouped_unscaled_demand.groupby("time.year").apply(
        _scale_demand, annual_demand=annual_demand_twh
    )

    # * 1e6 for TWh -> MWh
    return scaled_demand_profiles.to_array("end_use") * 1e6 * model_scaling_factor


def _scale_demand(
    one_year_profile: xr.Dataset, annual_demand: xr.Dataset
) -> xr.Dataset:
    "Scale demand in each year and sum different building types into one profile."
    year = one_year_profile.time.dt.year[0]
    normalised_profile = one_year_profile / one_year_profile.sum("time")
    demand = normalised_profile * annual_demand.sel(year=year).drop("year")
    return demand.sum("cat_name")


def prepare_annual_demand(annual_demand: pd.Series) -> xr.DataArray:
    """Restructure annual demand MultiIndex series into a multi-dimensional array.

    Result sums over all building categories and only contains hot water and space heating demands (not cooking).
    """
    return (
        annual_demand.rename_axis(columns="id")
        .stack()
        .unstack("end_use")
        .to_xarray()[["space_heat", "hot_water"]]
    )


def electrify_heat_demand_profiles(
    heat_demand: xr.DataArray, cop: xr.DataArray, electrification_shares: dict
) -> xr.DataArray:
    normalised_heat_demand = heat_demand / heat_demand.sum("time")

    # Weight COP based on _when_ it is demanded, A.K.A. the seasonal coefficient of performance (SCOP).
    weighted_average_cop = (
        (normalised_heat_demand * cop).sum("time").expand_dims(tech=["heat-pump"])
    )
    direct_electrification_eff = xr.DataArray([1], dims={"tech": ["direct"]})
    efficiency_da = xr.concat(
        [weighted_average_cop, direct_electrification_eff], dim="tech"
    )
    electrification_shares_da = (
        pd.DataFrame(electrification_shares)
        .rename_axis(index="tech")
        .to_xarray()
        .to_array("end_use")
    )
    assert np.isclose(
        electrification_shares_da.sum("tech"), 1
    ).all(), "Heat electrification shares must add up to 1."
    electrified_heat_demand = (
        heat_demand * electrification_shares_da / efficiency_da
    ).sum("tech")
    return electrified_heat_demand


if __name__ == "__main__":
    annual_demand = pd.read_csv(snakemake.input.annual_demand, index_col=[0, 1, 2])
    annual_demand_ds = prepare_annual_demand(annual_demand)
    unscaled_profiles = xr.open_dataset(snakemake.input.timeseries_data)
    scaled_profiles = scale_heat_demand_profiles(
        annual_demand_ds,
        unscaled_profiles,
        snakemake.params.sfh_mfh_shares,
        snakemake.params.scaling_factor,
    )
    if snakemake.wildcards.input_dataset == "electrified-heat":
        cop = xr.open_dataset(snakemake.input.cop).to_array("end_use")
        scaled_profiles = electrify_heat_demand_profiles(
            scaled_profiles, cop, snakemake.params.electrification_shares
        )

    final_df = (
        scaled_profiles.sum("end_use").astype("float32").to_series().unstack("id")
    )
    final_df.rename_axis(index="timesteps").to_csv(snakemake.output[0])
