# The workflow that builds models

When the pre-built models provide not enough flexibility for you, join us in the Euro-Calliope workshop where you can build your own models.
Next to pre-built models, Euro-Calliope contains the routines that allow to automatically build models from raw data.
All routines glued together form a scientific workflow.
When you execute the workflow yourself, you have more flexibility and more [customisation options](./customisation.md).
You should have experience executing scientific workflows or the willingness to learn.
The build process will require an internet connection, several gigabytes of available disk space, and a few hours of compute time.

## Prepare

!!! Info
    The workflow is developed on macOS, tested on macOS and Linux, but it cannot run natively on Windows.

1. Download your copy of the latest Euro-Calliope release: [![workflow DOI](https://img.shields.io/badge/workflow-10.5281/zenodo.3949793-blue)](https://doi.org/10.5281/zenodo.3949793) (or [clone any version from GitHub](https://github.com/calliope-project/euro-calliope)).

2. Create a conda environment from within which you can build the model.
This requires that you have [conda](https://conda.io) or [mamba](https://mamba.readthedocs.io/) installed on your machine.
We recommend `mamba`, as it's a faster drop-in replacement for `conda`.
Using either one, you can create the environment:

        # using mamba
        mamba env create -f environment.yaml --no-default-packages
        conda activate euro-calliope
        snakemake --list # test your installation

        # using conda
        conda env create -f environment.yaml --no-default-packages
        conda activate euro-calliope
        snakemake --profile profiles/conda --list # test your installation

3. Install a Gurobi license on your computer ([academic license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) comes at no cost), or [choose a different solver](./customisation.md#manual).

4. Create an account at the Copernicus Climate Data Service and a `$HOME/.cdsapirc` file with your credentials; see their [How To](https://cds.climate.copernicus.eu/api-how-to) (you do not need to manually install the client).

## Build

By default, the entire model will be built in the `./build/model` folder when you execute the workflow.
Execute the workflow like so:

```bash
snakemake
```

## Run

The build step creates all individual components of Euro-Calliope, like technologies and time series.
These can be combined to eventually build a final model to run simulations with.
For an example of such a model, see `./build/models/{resolution}/example-model.yaml`.
It is a complete Calliope model and can be used like any other, for example like so:

```bash
calliope run ./build/models/national/example-model.yaml
```

For more information on how to use Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.10/).

## Customise

Without any modifications, the workflow will rebuild the [pre-built models](../model/pre-built.md).
These models are examples and very likely require customisation to fit your purpose.
Once you've managed to build and run them, it's a good point in time to learn about [model customisation options](../model/customisation.md) and [workflow customisation options](./customisation.md).
