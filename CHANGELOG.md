# Changelog

## 1.1 (unpublished)

* **ADD** ability to move working directory (#45).
* **ADD** schema that automatically validates configuration files (#45).

* **UPDATE** JRC hydro database v4 -> v7 (#48). This entails two patches being removed:
    1. Romanian PHS data is no longer manually added from Geth et al. 2015.
    2. There are no index duplicates to drop in the JRC hydro database.

* **UPDATE** dependencies (#44):
    * Python 3.7 -> 3.8
    * atlite -> 0.2.1
    * geo packages from gdal 2.4 -> 3.2.1
    * Updates to NumPy, Pandas, xarray, and others

## 1.0.0 (2020-07-01)

First public release.
