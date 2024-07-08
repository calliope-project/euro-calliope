import math

import pandas as pd
import xarray as xr
from eurocalliopelib import utils


def scale_heat_demand_profiles(
    annual_demand: xr.Dataset,
    unscaled_demand_profiles: xr.Dataset,
    sfh_mfh_shares: dict,
) -> xr.DataArray:
    """Create a timeseries of demands for space heat and hot water demands across all building types.

    ASSUME: if calculating historic electrical heating requirements, annual electricity demand is distributed using the demand profile.
    It therefore ignores heat pump COP profiles and the possibility of a storage buffer.

    Args:
        annual_demand (xr.Dataset):
            Annual heat demand by year, resolution ID, end-use, and building category.
        unscaled_demand_profiles (xr.Dataset):
            Hourly heat demand profiles per year by resolution ID, end-use, and building category.
            Profiles will be normalised before scaling with annual demand, so their absolute magnitudes will be ignored.
        sfh_mfh_shares (dict):
            Share of single- and multi-family households, used to combine respective unscaled demand profiles into one "household" profile.
    Returns:
        xr.DataArray: `unscaled_demand_profiles` merged across building types and end uses, and scaled to have an annual sum equal to `annual_demand`.
    """
    assert math.isclose(
        sum(sfh_mfh_shares.values()), 1
    ), "Household type (single- vs multi-family home) ratios must add up to 1."

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
        _scale_demand, annual_demand=annual_demand
    )

    return scaled_demand_profiles.to_array("end_use").sum("end_use")


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


if __name__ == "__main__":
    annual_demand = pd.read_csv(snakemake.input.annual_demand, index_col=[0, 1, 2])
    annual_demand_ds = prepare_annual_demand(annual_demand)
    unscaled_profiles = xr.open_dataset(snakemake.input.timeseries_data)
    scaled_profiles = scale_heat_demand_profiles(
        annual_demand_ds, unscaled_profiles, snakemake.params.sfh_mfh_shares
    )

    # Demands are stored as negative values for Calliope to ingest
    if "historic" not in snakemake.wildcards.input_dataset:
        scaled_profiles *= -1

    final_df = scaled_profiles.astype("float32").to_series().unstack("id")
    final_df.to_csv(snakemake.output[0])
