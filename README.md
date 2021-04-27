# Euro Calliope

A model of the European electricity system built using Calliope.

This repository contains the workflow routines that automatically build models from source data. Alternatively to building models yourself, you can use [pre-built models](https://doi.org/10.5281/zenodo.3949553) that run out-of-the-box. You can find a more detailed description of the first application in a [scientific article in Joule](https://doi.org/10.1016/j.joule.2020.07.018).

[![article DOI](https://img.shields.io/badge/article-10.1016/j.joule.2020.07.018-blue)](https://doi.org/10.1016/j.joule.2020.07.018)
[![pre-built models DOI](https://img.shields.io/badge/prebuilts-10.5281%2Fzenodo.3949553-blue)](https://doi.org/10.5281/zenodo.3949553)

## At a glance

euro-calliope models the European electricity system with each location representing an administrative unit. It can be built on three spatial resolutions: on the continental level as a single location, on the national level with 34 locations, and on the regional level with 497 locations. On each node, renewable generation capacities (wind, solar, bioenergy) and balancing capacities (battery, hydrogen) can be built. In addition, hydro electricity and pumped hydro storage capacities can be built up to the extent to which they exist today. All capacities are used to satisfy electricity demand on all locations which is based on historic data. Locations are connected through transmission lines of unrestricted capacity. Using [Calliope](https://www.callio.pe), the model is formulated as a linear optimisation problem with total monetary cost of all capacities as the minimisation objective. All elements of euro-calliope can be manipulated either by changing the configuration in `config/default.yaml` or by manipulating the build workflow before building the model.

## Get ready to build the model

1. The workflow is developed on macOS, should run on Linux (untested), and cannot run natively on Windows.

2. You need [conda](https://conda.io/docs/index.html) to build and use the model. Using conda, you can create a conda environment from within you can build the model:

```
conda env create -f environment.yaml
conda activate euro-calliope
```

3. You need a Gurobi license installed on your computer, or you need to choose a different solver.

4. You need an account at the Copernicus Climate Data Service and you need to create a `$HOME/.cdsapirc` file with your credentials, see their [How To](https://cds.climate.copernicus.eu/api-how-to) (you do not need to manually install the client).

5. Further, you need all data files that cannot be retrieved automatically:

* [Maritime Boundaries v10 -> World Exclusive Economic Zones v10](http://www.marineregions.org/downloads.php), to be placed in `./data/World_EEZ_v10_20180221`

## Build the model

Because input data is large, the actual model including its data is not part of this repository. To use the model, you first need to build it from input data and scripts. Running the build step will build the model in the `./model` folder.

    snakemake --use-conda

## Build the model on a cluster

You may want to build the model on a cluster. While you can build euro-calliope on [any cluster that is supported by Snakemake](https://snakemake.readthedocs.io/en/stable/executing/cluster.html), our default configuration is targeted at ETH's Euler cluster. To build the model on Euler, use the following command:

    snakemake --use-conda --profile config/euler

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Work local, build on remote

If you are like us, you may want to work locally (to change configuration parameters, add modules etc), but execute remotely on the cluster. We support this workflow through three Snakemake rules: `send`, `receive`, and `clean_cluster_results`. It works like the following.

First, start local and make sure the `cluster-sync` configuration parameters fit your environment. Next, run `snakemake --use-conda send` to send the entire repository to your cluster. On the cluster, execute the workflow with Snakemake (see above). After the workflow has finished, download results by locally running `snakemake --use-conda receive`. By default, this will download results into `build/euler`.

This workflow works iteratively too. After analysing your cluster results locally, you may want to make changes locally, send these changes to the cluster (`snakemake --use-conda send`), rerun on the cluster, and download updated results (`snakemake --use-conda receive`).

To remove cluster results on your local machine, run `snakemake --use-conda clean_cluster_results`.

## Example use of the model

The build step creates all individual components of `euro-calliope`, like technologies and time series. These can be combined to eventually build a final model to run simulations with. For an example of such a model, see `./build/model/{resolution}/example-model.yaml`. It is a complete Calliope model and can be used like any other, for example like this:

```Bash
$ calliope run ./build/model/national/example-model.yaml
```

For more information on how to use and modify Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io).

## Model components

After a successful full build (see "Build the model"), the following files will exist in `build/model`:

```
├── {resolution}                           <- For each spatial resolution an individual folder.
│   ├── capacityfactors-{technology}.csv   <- Timeseries of capacityfactors of all renewables.
│   ├── directional-rooftop.yaml           <- Override discriminating rooftop PV by orientation.
│   ├── electricity-demand.csv             <- Timeseries of electricity demand on each node.
│   ├── example-model.yaml                 <- Calliope model definition.
│   ├── link-all-neighbours.yaml           <- Connects neighbouring locations with transmission.
│   ├── locations.csv                      <- Map from Calliope location id to name of location.
│   └── locations.yaml                     <- Defines all locations and their max capacities.
├── build-metadata.yaml                    <- Metadata of the build process.
├── demand-techs.yaml                      <- Definition of demand technologies.
├── environment.yaml                       <- A conda file defining an environment to run the model.
├── interest-rate.yaml                     <- Interest rate of all capacities.
├── link-techs.yaml                        <- Definition of link technologies.
├── README.md                              <- Documentation.
├── renewable-techs.yaml                   <- Definition of supply technologies.
└── storage-techs.yaml                     <- Definition of storage technologies.
```

Alternatively to a full build, each of these model components can be built individually, by running `snakemake --use-conda <path-to-component>`. The model components can be used to [configure a Calliope model](https://calliope.readthedocs.io/en/stable/user/building.html). For an example model configuration, see "Example use of the model" above.

## Units and scaling

The default units used within euro-calliope are `100 GW`, `100 GWh`, `billion EUR`, and `10,000 km2`. All data going into Calliope and all Calliope result data will be given using these units. While they may be unusual, these units lead to a numerical model that is well suited for the interior-point solution algorithm that is used by default. The units are tuned so as to work best for models with a time resolution of a few hours and a duration of one year. For other types of problems, or other solution algorithms, the units may need to be changed to avoid numerical issues within the solver.

You can easily change the units and scale all values using the `scaling-factor` configuration values in `config/default.yaml`. However, these values must be changed before building the model. You may want to run `snakemake clean` before changing these values. The base units on which the scaling factors are applied are `1 MW`, `1 MWh`, `EUR`, and `km2`. So for example, the default unit for energy (100 GWh) is derived by scaling the base unit (1 MWh) with a scaling factor of `0.00001`.

## Repo structure

* `build/model`: Contains the entire model after the build step, including Calliope definition files and data (does not exist initially).
* `config`: Files with configuration parameters which influence the model build process.
* `data`: Small input data used within the model build process.
* `docs`: Documentation of the model and the build process.
* `envs`: Files defining the conda environment which are used to build the model.
* `lib`: Library code in form of a Python package that is reused in many places of this repository.
* `notebooks`: Notebooks for various data analysis or preparation steps. Not within main workflow.
* `rules`: Snakemake workflows defining the build process.
* `scripts`: Contains the scripts to build the model.
* `templates`: Contains templates of Calliope model files.
* `tests`: Contains a test usage of the model.

## Run the tests

Tests of models with continental and national resolution run automatically when you run the entire workflow. To run the tests of models with regional resolution do the following:

    snakemake --use-conda build/logs/regional/test-report.html

Exchanging `regional` with `national` or `continental` allows you to run tests on the respective resolution explicitly.

## Run minimal test

As a developer, you may want to run the entire workflow often to spot errors early. For that, you can use a minimal test configuration that takes less time to run.

    snakemake --use-conda --configfile="config/minimal.yaml"

Make sure to run this in a clean working directory. Do not use the working directory in which you are using your normal configuration.

## License

euro-calliope is developed and maintained within the [Calliope project](https://www.callio.pe). The code in this repository is MIT licensed.
