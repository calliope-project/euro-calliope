from functools import partial
from multiprocessing import Pool

import pandas as pd
import xarray as xr


def group_gridcells(
    gridded_data: xr.Dataset,
    grid_weight: xr.DataArray,
    annual_demand: xr.DataArray,
    threads: int,
) -> pd.DataFrame:
    """Group gridded heat data into resolution-specific units, taking a weighted average of grid values.

    Also take a weighted average of hot water and space heat demand (using per-unit annual demands) to return a single "heat" timeseries.

    Args:
        gridded_data (xr.Dataset): Gridded timeseries space heat and hot water data.
        grid_weight (xr.DataArray): Weighted mapping from grid (a.k.a. "site") to units.
        annual_demand (xr.DataArray): Per-unit annual space heating and hot water demands.
        threads (int): Number of threads over which to undertake multiprocessing.

    Returns:
        pd.DataFrame: Index: timeseries, Columns: unit IDs
    """

    partial_func = partial(
        _site_weighted_ave, gridded_data=gridded_data, grid_weight=grid_weight
    )
    # This is a slow operation, so we parallelise it.
    with Pool(threads) as pool:
        per_id_averages = pool.map(partial_func, grid_weight.id.values)
    weighted_average_ds = xr.concat(per_id_averages, dim="id")

    weighted_average_da = (
        weighted_average_ds.to_array("end_use")
        .groupby("time.year")
        .apply(_end_use_weighted_ave, annual_demand=annual_demand)
    )
    return weighted_average_da.astype("float32").to_series().unstack("id")


def prepare_annual_demand(annual_demand: pd.Series) -> xr.DataArray:
    """Restructure annual demand MultiIndex series into a multi-dimensional array.

    Result sums over all building categories and only contains hot water and space heating demands (not cooking).
    """
    return (
        annual_demand.rename_axis(columns="id")
        .stack()
        .to_xarray()
        .sel(end_use=["space_heat", "hot_water"])
        .sum("cat_name")
    )


def _site_weighted_ave(
    id: str, gridded_data: xr.Dataset, grid_weight: xr.DataArray
) -> xr.Dataset:
    """Multi-processing helper function to get the weighted average of all gridcells for a given spatial unit (id)."""
    id_grid_weight = grid_weight.sel(id=id).dropna("site")
    return (
        gridded_data.reindex_like(id_grid_weight)
        # we normalise the weights in case they aren't already
        * (id_grid_weight / id_grid_weight.sum("site"))
    ).sum(["site"])


def _end_use_weighted_ave(
    one_year_profile: xr.DataArray, annual_demand: xr.DataArray
) -> xr.DataArray:
    """Take a weighted average of all heat energy end uses using annual demands per spatial unit as the weights."""
    year = one_year_profile.time.dt.year[0]
    normalised_demand = annual_demand / annual_demand.sum("end_use")
    demand = one_year_profile * normalised_demand.sel(year=year).drop("year")
    return demand.sum("end_use")


if __name__ == "__main__":
    gridded_data = xr.open_dataset(snakemake.input.gridded_timeseries_data)
    grid_weights = xr.open_dataarray(snakemake.input.grid_weights)
    annual_demand = pd.read_csv(snakemake.input.annual_demand, index_col=[0, 1, 2])
    annual_demand_ds = prepare_annual_demand(annual_demand)
    resolution_specific_data = group_gridcells(
        gridded_data, grid_weights, annual_demand_ds, snakemake.threads
    )
    resolution_specific_data.to_csv(snakemake.output[0])
