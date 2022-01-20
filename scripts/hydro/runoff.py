"""Retrieve runoff data from ERA5 as atlite cutout."""
import logging

import atlite


def runoff(path_to_cutout, first_year, final_year, x_min, x_max, y_min, y_max):
    """Retrieve runoff data from ERA5 as atlite cutout."""
    logging.basicConfig(level=logging.INFO)
    x_range = slice(x_min, x_max)
    y_range = slice(y_max, y_min)
    time_range = slice(f"{first_year}-01", f"{final_year}-12")

    cutout = atlite.Cutout(
        path=path_to_cutout,
        module="era5",
        xs=x_range,
        ys=y_range,
        time=time_range
    )
    cutout.prepare(["runoff"])


if __name__ == "__main__":
    runoff(
        first_year=snakemake.params.first_year,
        final_year=snakemake.params.final_year,
        x_min=snakemake.params.x_min,
        x_max=snakemake.params.x_max,
        y_min=snakemake.params.y_min,
        y_max=snakemake.params.y_max,
        path_to_cutout=snakemake.output[0]
    )
