"""Retrieve runoff data from ERA5 as atlite cutout."""
from pathlib import Path

import atlite
import xarray as xr

X_RANGE = slice(-10.7, 34.6)
Y_RANGE = slice(71.2, 34.5)


def runoff(path_to_cutout, year):
    """Retrieve runoff data from ERA5 as atlite cutout."""
    path_to_cutout = Path(path_to_cutout)

    monkeypatch_atlite_to_download_runoff_only()

    cutout = atlite.Cutout(
        name=path_to_cutout.name,
        module="era5",
        cutout_dir=path_to_cutout.parent,
        xs=X_RANGE,
        ys=Y_RANGE,
        years=slice(year - 1, year),
        months=slice(1, 12)
    )
    cutout.prepare()


def prepare_month_era5(year, month, xs, ys):
    area = atlite.datasets.era5._area(xs, ys)

    # Reference of the quantities
    # https://confluence.ecmwf.int/display/CKB/ERA5+data+documentation
    # (shortName) | (name)                                      | (paramId)
    # tisr        | TOA incident solar radiation                | 212
    # ssrd        | Surface Solar Rad Downwards                 | 169
    # ssr         | Surface net Solar Radiation                 | 176
    # fdir        | Total sky direct solar radiation at surface | 228021
    # ro          | Runoff                                      | 205
    # 2t          | 2 metre temperature                         | 167
    # sp          | Surface pressure                            | 134
    # stl4        | Soil temperature level 4                    | 236
    # fsr         | Forecast surface roughnes                   | 244

    with atlite.datasets.era5._get_data(area=area, year=year, month=month, variable=['runoff']) as ds, \
            atlite.datasets.era5._get_data(area=area, year=year, month=month, day=1, variable=['orography']) as ds_m:
        ds_m = ds_m.isel(time=0, drop=True)
        ds = xr.merge([ds, ds_m], join='left')

        ds = atlite.datasets.era5._rename_and_clean_coords(ds)
        ds = atlite.datasets.era5._add_height(ds)

        ds = ds.rename({'ro': 'runoff', })

        yield (year, month), ds


def monkeypatch_atlite_to_download_runoff_only():
    # this saves ~90% of data to download and store
    atlite.datasets.era5.weather_data_config = {
        '_': dict(
            tasks_func=atlite.datasets.era5.tasks_monthly_era5,
            prepare_func=prepare_month_era5
        )
    }


if __name__ == "__main__":
    runoff(
        year=snakemake.params.year,
        path_to_cutout=snakemake.output[0]
    )
