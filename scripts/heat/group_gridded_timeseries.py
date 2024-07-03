from functools import partial
from multiprocessing import Pool

import pandas as pd
import xarray as xr


def group_gridcells(
    gridded_data: xr.Dataset, grid_weight: xr.DataArray, threads: int
) -> xr.Dataset:
    """Group gridded heat data into resolution-specific units, taking a weighted average of grid values.

    Args:
        gridded_data (xr.Dataset): Gridded timeseries space heat and hot water data.
        grid_weight (xr.DataArray): Weighted mapping from grid (a.k.a. "site") to units.
        threads (int): Number of threads over which to undertake multiprocessing.

    Returns:
        xr.Dataset: `gridded_data` with `sites` dimension reduced to an `id` dimension, representing resolution-specific units.
    """

    apply_weights = partial(
        _site_weighted_ave, gridded_data=gridded_data, grid_weight=grid_weight
    )
    # This is a slow operation, so we parallelise it.
    with Pool(threads) as pool:
        per_id_averages = pool.map(apply_weights, grid_weight.id.values)
    weighted_average_ds = xr.concat(per_id_averages, dim="id")

    return weighted_average_ds


def group_end_uses(
    resolution_specific_data: xr.Dataset, annual_demand: xr.DataArray
) -> xr.DataArray:
    """Take a weighted average of hot water and space heat end uses (using per-unit annual demands) to return a single "heat" end use timeseries.

    Args:
        resolution_specific_data (xr.Dataset): Timeseries space heat and hot water data that has been pre-grouped into resolution-specific units.
        annual_demand (xr.DataArray): Per-spatial unit annual space heating and hot water demands.

    Returns:
        xr.DataArray: `resolution_specific_data` with end uses averaged into a single `heat` end use.
    """
    weighted_average_da = (
        resolution_specific_data.to_array("end_use")
        .groupby("time.year")
        .apply(_end_use_weighted_ave, annual_demand=annual_demand)
    )
    return weighted_average_da


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
    """Get the weighted average of all gridcells for a given spatial unit (id).

    This function exists to enable multi-processing across IDs.
    """
    id_grid_weight = grid_weight.sel(id=id).dropna("site")
    normalised_weight = id_grid_weight / id_grid_weight.sum("site")
    return (gridded_data.reindex_like(id_grid_weight) * normalised_weight).sum("site")


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
        gridded_data, grid_weights, snakemake.threads
    )
    resolution_specific_grouped_end_use_data = group_end_uses(
        resolution_specific_data, annual_demand_ds
    )
    final_df = (
        resolution_specific_grouped_end_use_data.astype("float32")
        .to_series()
        .unstack("id")
    )
    final_df.to_csv(snakemake.output[0])
