# Release Notes

## 2.0 (unpublished)

* **ADD** Pipeline to download and process JRC-IDEES datasets for use in industry/tertiary/transport sector rules (#133)
* **ADD** Download and preliminary pre-processing of Eurostat and equivalent Swiss data, for use in most downstream sector-coupling processing (#113)

## 1.1 (unpublished)

* **ADD** documentation to ReadTheDocs (#114) and extend it massively:
    * **ADD** an entire section for each use case: the pre-built models and the workflow (#154).
    * **ADD** a troubleshooting section (#159).
    * **ADD** templates to our issue tracker (#156).
    * **ADD** a visualisation of directed acyclic graph spanned by all workflow rules (#108).
* **ADD** Link YAML file defining national net-transfer capacities from ENTSO-E TYNDP 2020 scenarios (#61).
* **ADD** automatic download of EEZ (#99).
* **ADD** Danish Energy Agency and Schroeder et al (2013) cost data as well as no hydro fixed costs as optional overrides (#18, #129).
* **ADD** override to constrain energy to power ratios of battery and hydrogen storage (#130).
* **ADD** ability to move working directory (#45).
* **ADD** schema that automatically validates configuration files (#45).
* **ADD** automatic download of hydro basins data (#34).
* **ADD** minimal configuration to be able to test the entire workflow more quickly (#60).
* **ADD** installation of `curl` and `unzip` from conda-forge, to increase portability (#59).
* **ADD** sync infrastructure to easily send and receive files to and from a cluster (#74).
* **ADD** parameter `station-nearest-basin-max-km` controlling the mapping of hydro power stations to basins (#138).
* **ADD** optional email notifications whenever builds fail or succeed (#92).
* **ADD** option to shed load at all locations (#131).
* **ADD** option to choose solar and wind potential scenario, limiting the amount of eligible surfaces (#153).

* **UPDATE** Calliope version from 0.6.5 to [0.6.7](https://calliope.readthedocs.io/en/v0.6.7/history.html#id1) (#73).
* **UPDATE** EEZ updated from v10 to v11 (difference in offshore area is < 1% for all relevant countries) (#99).
* **UPDATE** ENTSO-E national electricity load data gap filling methods (priority order, interpolation distance, outlier handling, 29th Feb handling) included in config (#42, #91).
* **UPDATE** IRENA hydro generation data from 2018 to 2020 (#40).
* **UPDATE** Make scaling pumped hydro capacity according to Geth et al. (2015) optional with a boolean config parameter `scale-phs-according-to-geth-et-al` which defaults to False (no scaling) (#49).
* **UPDATE** JRC hydro database v4 -> v9 (#48, #57). This entails one patch being removed:
    1. Romanian PHS data is no longer manually added from Geth et al. (2015).
* **UPDATE** Improve gap-filling method for national electricity load data (#3).
* **UPDATE** dependencies (#44, #73, #142):
    * Python 3.7 -> 3.8
    * Snakemake 5.8.2 -> 6.1.1 (uses mamba by default)
    * atlite -> 0.2.1
    * geo packages from gdal 2.4 -> 3.2.1
    * Updates to NumPy, Pandas, xarray, pytest, and others

* **FIX** the centroid determination of all locations which had been calculated on an unprojected reference system before and was therefore slightly off (#147).

## 1.0.0 (2020-07-01)

First public release.
