# Euro Calliope

A model of the European power system built using Calliope.

This repository contains the model and the routines that automatically generate the model from source data.

## Get ready

You need [conda](https://conda.io/docs/index.html) to use the model. Using conda, you can create a conda environment from within you can run it:

    conda env create -f requirements.yml

Don't forget to activate the environment.

## Example use of the model

Coming soon.

## Repo structure

* `model`: contains the entire model, including Calliope definition files and data
* `src`: contains the source data and source code to generate the model
* `tests`: contains a test usage of the model

## Generate the model

The following applies to you if you want to manipulate the structure of Euro Calliope itself. It does not apply to you when you intend to _use_ Euro Calliope.

### Get ready

You need [conda](https://conda.io/docs/index.html) to generate the model. Using conda, you can create a conda environment from within you can generate it:

    conda env create -f src/envs/default.yaml

Don't forget to activate the environment.

Further, you need all data files that cannot be retrieved automatically:

* country shapes and their national electricity demand, to be placed in `./src/data/national-technical-potential.geojson` # FIXME should come from Zenodo
* shapes of exclusive economic zones (eez), to be placed in `./src/data/eez-in-europe.geojson` # FIXME should come from Zenodo
* fraction of shared coasts, nations to eez, to be placed in `./src/data/national-shared-coast.csv`
* national land eligibility of renewables, to be placed in `./src/data/national-eligibility.csv` # FIXME should come from Zenodo
* spatio-temporal capacity factors in `./src/data/capacityfactors/`, where time and space dimensions are defined by two files: # FIXME should come from Zenodo
    * an id map, where each pixel points to a time series: `./src/data/capacityfactors/{technology}-ids.tif`
    * all indexed time series: `./src/data/capacityfactors/{technology}-timeseries.nc`

### Regenerate the model

Running the generation step will override files in the `./model` folder. Make sure this is what you want to do and that you have a copy.

    snakemake --use-conda

### Run the tests

    snakemake test --use-conda
