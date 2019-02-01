# Euro Calliope

A model of the European power system built using Calliope.

This repository contains the routines that automatically generate the model from source data.

## Get ready

1. You need [conda](https://conda.io/docs/index.html) to build and use the model. Using conda, you can create a conda environment from within you can run it:

    conda env create -f requirements.yml
    conda activate euro-calliope

2. Because of many recent changes in the way conda works, there is currently a [compatibility issue](https://bitbucket.org/snakemake/snakemake/issues/1029/subshells-and-conda-44) between `Snakemake 5.4` and recent conda versions (at least starting from `conda 4.6`). A workaround is to have the `<path-to-conda-root>/bin` folder on your `PATH` (make sure the folder contains `activate`). Because this is discouraged by the conda devs, it is a good idea to add it to `PATH` only when calling `snakemake` -- for example in a function ([my setup as example](https://github.com/timtroendle/.settings/blob/a5afc0c5f37afe4f5b1b924639e03c130fc7bdb7/fish/functions/smake.fish#L1)).

3. Further, you need all data files that cannot be retrieved automatically:

* country shapes and their national electricity demand, to be placed in `./src/data/national-technical-potential.geojson` # FIXME should come from Zenodo
* shapes of exclusive economic zones (eez), to be placed in `./src/data/eez-in-europe.geojson` # FIXME should come from Zenodo
* fraction of shared coasts, nations to eez, to be placed in `./src/data/national-shared-coast.csv`
* national land eligibility of renewables, to be placed in `./src/data/national-eligibility.csv` # FIXME should come from Zenodo
* spatio-temporal capacity factors in `./src/data/capacityfactors/`, where time and space dimensions are defined by two files: # FIXME should come from Zenodo
    * an id map, where each pixel points to a time series: `./src/data/capacityfactors/{technology}-ids.tif`
    * all indexed time series: `./src/data/capacityfactors/{technology}-timeseries.nc`

## Generate the model

Because input data is large, the actual model including this data is not part of this repository. To use the model, you first need to generate it from input data and scripts. Running the generation step will generate the model in the `./model` folder.

    snakemake --use-conda

## Run on Euler cluster

To run on Euler, use the following command:

    snakemake --use-conda --profile config/euler

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/stable/executable.html#cluster-execution) and take `config/euler` as a starting point.

## Example use of the model

The generation step created all single parts of `euro-calliope`, like technologies and time series. These can be combined to eventually build a final model to run simulations with. For an example of such a model, see `./tests/simple-model.yaml`. It is a complete Calliope model and can be used like any other, for example like this:

```Bash
$ calliope run ./tests/simple-model.yaml
```

For more information on how to use Calliope models, see [Calliope's documentation](https://www.callio.pe).

## Repo structure

* `model`: contains the entire model after the generation step, including Calliope definition files and data
* `src`: contains the source data and source code to generate the model
* `tests`: contains a test usage of the model

## Run the tests

    snakemake test --use-conda
