# Release Notes

## sector-coupled (unpublished)

### Added (workflow)
* **ADD** Pipeline to download and process JRC-IDEES datasets for use in industry/tertiary/transport sector rules (#133)
* **ADD** Download and preliminary pre-processing of Eurostat and equivalent Swiss data, for use in most downstream sector-coupling processing (#113)

## 1.2.0 (unpublished)

### Added (models)

* **ADD** top level model configuration file `scenarios.yaml` which can act as the focal point to consolidate all available overrides in the model (#209).
* **ADD** nuclear power plant technology with capacity limits. Capacity limits can be equals to today or be bound by a minimum and maximum capacity to represent an available range in future. In either case, capacities are allocated at a subnational resolution based on linear scaling from current capacity geolocations, using the JRC power plant database (#78).

### Added (workflow)

* **ADD** configuration option to build model timeseries data over multiple years, using `first-year` and `final-year` temporal scopes. Available years are 2010-2016 at time of implementing functionality (#152).
* **ADD** nuclear technology capacity allocation workflow which uses the configuration parameter `nuclear-capacity-scenario` to select whether today's capacities define limits in the model definition ("current") or whether ranges set bounds on future capacity (by linking to a configuration CSV file) (#78).
### Updated (models)

* **UPDATED** Final model configuration and data files structure (#145) to:
    * make each spatial resolution model self-contained (i.e., no shared files between resolutions);
    * split technology definitions into self-explanatory files and into subdirectories named after Calliope abstract technology groups (e.g., `supply/wind-offshore.yaml` for offshore wind supply technology). This enables technologies to be added to or removed from the model by simply changing the model configuration file import list.;
    * keep technology definitions and their allocations to locations in the model in the same file; and
    * separate tech config YAML files from data CSV files. The former are found in the `techs` subdirectory, while the latter are in `timeseries`.

### Updated (workflow)

* **UPDATED** YAML templates and parametrisation restructured:
    * Parametrisation moved to eurocalliopelib.
    * Rules to parametrise split into smaller technology-specific rules, to ensure inputs are directly relevant to the files being parametrised.
    * YAML templates restructured to match structure of final model (see `Updated (models) above`);

* **UPDATE** cluster sync infrastructure to retain file permission defaults on the cluster. This change improves team collaboration, as default group settings will apply to the files on the cluster (#214).

### Fixed (models)

* **FIX** spatial proxy of `landscape-care-residues` biofuel in the default configuration. National potentials were spatially allocated to sub-national regions based on `population` rather than `forest` land use. As the potential of `landscape-care-residues` is low, this change has only a minor impact on the model. Continentally, only 3.5% of the total biofuel potential is of this type. Nationally, only three countries have a share of slightly above 10%: Bosnia and Herzegovina, Albania, and Norway. National potentials are unaffected by this change.

### Fixed (documentation)

* **FIX** links in the documention to always point to the most recent version of the pre-builts (#218).

## 1.1.0 (2021-12-22)

This version is a minor update to v1.0. It comes with updated input data, updated dependencies, several convenience features, several optional overrides, and with massively extended documentation. You should not expect model results to deviate much from v1.0.

In the following, we split additions, fixes, and updates between those which affect the data and configuration in the final _models_, and those which affect the _workflow_ creating the models. If you are a user of the pre-built models, you only need to care about the former. If you build models yourself, you likely care about all changes.

### Added (models)

* **ADD** documentation to ReadTheDocs (#114) and extend it massively:
    * **ADD** an entire section for each use case: the pre-built models and the workflow (#154).
    * **ADD** a troubleshooting section (#159).
    * **ADD** templates to our issue tracker (#156).
    * **ADD** a visualisation of directed acyclic graph spanned by all workflow rules (#108).
    * **ADD** a complete enumeration of the configuration parameters of the workflow (#68).
    * **ADD** a developer documentation that aids contributions (#109).
    * **ADD** a visualisation of administrative units on all three spatial resolutions (#165).
* **ADD** national net-transfer capacities from ENTSO-E TYNDP 2020 scenarios as optional override (#61).
* **ADD** Danish Energy Agency and Schroeder et al (2013) cost data as well as no hydro fixed costs as optional overrides (#18, #129).
* **ADD** override to constrain energy to power ratios of battery and hydrogen storage (#130).
* **ADD** option to shed load at all locations (#131).

### Added (workflow)

* **ADD** automatic downloads of the following datasources:
    * EEZ (#99),
    * hydro basins data (#34),
    * biofuels potential and cost data (#194).
* **ADD** ability to move working directory (#45).
* **ADD** schema that automatically validates configuration files (#45).
* **ADD** minimal configuration to be able to test the entire workflow more quickly (#60).
* **ADD** installation of `curl` and `unzip` from conda-forge, to increase portability (#59).
* **ADD** sync infrastructure to easily send and receive files to and from a cluster (#74).
* **ADD** parameter `station-nearest-basin-max-km` controlling the mapping of hydro power stations to basins (#138).
* **ADD** optional email notifications whenever builds fail or succeed (#92).
* **ADD** option to choose solar and wind potential scenario, limiting the amount of eligible surfaces (#153).

### Updated (models)

* **UPDATE** Calliope version from 0.6.5 to [0.6.7](https://calliope.readthedocs.io/en/v0.6.7/history.html#id1) (#73).
* **UPDATE** EEZ updated from v10 to v11 (difference in offshore area is < 1% for all relevant countries) (#99).
* **UPDATE** IRENA hydro generation data from 2018 to 2020 (#40).
* **UPDATE** JRC hydro database v4 -> v9 (#48, #57). This entails one patch being removed:
    1. Romanian PHS data is no longer manually added from Geth et al. (2015).
* **UPDATE** By default, pumped hydro capacity is based on JRC hydro database only, not according to Geth et al. (2015) (#49).
* **UPDATE** JRC biofuel potentials data source from @RuizCastello:2015 to ENSPRESO (@Ruiz:2019) (#194). This changes the potentials of the continent (0.1% less), Montenegro (93% less), North Macedonia (44% less), the UK (6% more), and the Netherlands (1% less).
* **UPDATE** Improve gap-filling method for national electricity load data (#3).

### Updated (workflow)

* **UPDATE** Workflow build metadata is built as part of the rule `all` instead of the rule `model` (#126).
* **UPDATE** ENTSO-E national electricity load data gap filling methods (priority order, interpolation distance, outlier handling, 29th Feb handling) included in config (#42, #91).
* **UPDATE** location of non-automatically derivable datasets: they are not included in the workflow repository anymore, but instead published separately on Zenodo (#201).
* **UPDATE** Make scaling pumped hydro capacity according to Geth et al. (2015) optional with a boolean config parameter `scale-phs-according-to-geth-et-al` which defaults to False (no scaling) (#49).
* **UPDATE** dependencies (#44, #73, #142):
    * Python 3.7 -> 3.8
    * Snakemake 5.8.2 -> 6.1.1 (uses mamba by default)
    * atlite -> 0.2.1
    * geo packages from gdal 2.4 -> 3.2.1
    * Updates to NumPy, Pandas, xarray, pytest, and others

### Fixed (models)

* **FIX** variable cost of biofuels. This reduces variable cost of biofuels from 64.83 to 44.14 €/MWh_el (32%) using the default settings.
    * Fix cost of municipal waste (#193).
    * Fix alignment of cost and potential years for biofuels (#195).
* **FIX** the centroid determination of all locations which had been calculated on an unprojected reference system before and was therefore slightly off (#147).

## 1.0.0 (2020-07-01)

First public release.
