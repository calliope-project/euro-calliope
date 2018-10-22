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

    conda env create -f src/envs/default.yml

Don't forget to activate the environment.

Further, you need all data files that cannot be retrieved automatically:

    * national technical potential of renewables, to be place in `./src/data/national-technical-potential.geojson`

### Regenerate the model

Running the generation step will override files in the `./model` folder. Make sure this is what you want to do and that you have a copy.

    snakemake --use-conda

### Run the tests

    snakemake test
