# Changelog

## 1.1 (unpublished)

* **ADD** ability to move working directory (#45).
* **ADD** schema that automatically validates configuration files (#45).
* **ADD** automatic download of hydro basins data (#34).
* **ADD** minimal configuration to be able to test the entire workflow more quickly (#60).
* **ADD** installation of `curl` and `unzip` from conda-forge, to increase portability (#59).
* **ADD** sync infrastructure to easily send and receive files to and from a cluster (#74).

* **UPDATE** ENTSOE national electricity load data gap filling priority order included in config (#42)
* **UPDATE** IRENA hydro generation data from 2018 to 2020 (#40).
* **UPDATE** Make scaling pumped hydro capacity according to Geth et al. (2015) optional with a boolean config parameter `scale-phs-according-to-geth-et-al` which defaults to False (no scaling) (#49).
* **UPDATE** JRC hydro database v4 -> v7 (#48). This entails one patch being removed:
    1. Romanian PHS data is no longer manually added from Geth et al. (2015).
* **UPDATE** Improve gap-filling method for national electricity load data (#3).
* **UPDATE** dependencies (#44):
    * Python 3.7 -> 3.8
    * atlite -> 0.2.1
    * geo packages from gdal 2.4 -> 3.2.1
    * Updates to NumPy, Pandas, xarray, and others

## 1.0.0 (2020-07-01)

First public release.
