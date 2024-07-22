from functools import partial
from multiprocessing import Pool

import xarray as xr


def group_gridcells(
    gridded_data: xr.Dataset, grid_weight: xr.DataArray, threads: int
) -> xr.DataArray:
    """Group gridded heat data into resolution-specific units, taking a weighted average of grid values.

    Args:
        gridded_data (xr.Dataset): Gridded timeseries space heat and hot water data.
        grid_weight (xr.DataArray): Weighted mapping from grid (a.k.a. "site") to units.
        threads (int): Number of threads over which to undertake multiprocessing.

    Returns:
        xr.DataArray: `gridded_data` with `sites` dimension reduced to an `id` dimension, representing resolution-specific units.
    """

    apply_weights = partial(
        _site_weighted_ave, gridded_data=gridded_data, grid_weight=grid_weight
    )
    # This is a slow operation, so we parallelise it.
    with Pool(threads) as pool:
        per_id_averages = pool.map(apply_weights, grid_weight.id.values)
    weighted_average_ds = xr.concat(per_id_averages, dim="id")

    return weighted_average_ds


def _site_weighted_ave(
    id: str, gridded_data: xr.Dataset, grid_weight: xr.DataArray
) -> xr.Dataset:
    """Get the weighted average of all gridcells for a given spatial unit (id).

    This function exists to enable multi-processing across IDs.
    """
    id_grid_weight = grid_weight.sel(id=id).dropna("site")
    normalised_weight = id_grid_weight / id_grid_weight.sum("site")
    return (gridded_data.reindex_like(id_grid_weight) * normalised_weight).sum("site")


if __name__ == "__main__":
    gridded_data = xr.open_dataset(snakemake.input.gridded_timeseries_data)
    grid_weights = xr.open_dataarray(snakemake.input.grid_weights)

    resolution_specific_data = group_gridcells(
        gridded_data, grid_weights, snakemake.threads
    )
    resolution_specific_data.to_netcdf(snakemake.output[0])
