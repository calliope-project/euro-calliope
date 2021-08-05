# Build the models yourself using the workflow

When the pre-built models provide not enough flexibility for you, join us in the euro-calliope workshop where you can build your own models. Next to pre-built models, Euro Calliope contains the routines that allow to automatically build models from raw data. All routines glued together form a scientific workflow. When you execute the workflow yourself, you have more flexibility and more [customisation options](./customisation). You should have experience executing scientific workflows or the willingness to learn. The build process will require an internet connection, several gigabytes of available disk space, and a few hours of compute time.

## Prepare

!!! note
    The workflow is developed on macOS, tested on macOS and Linux, but it cannot run natively on Windows.

1. [Download](https://doi.org/10.5281/zenodo.3949793) your copy of the latest euro-calliope release or [clone any version from GitHub](https://github.com/calliope-project/euro-calliope).

2. Create a conda environment from within which you can build the model. This requires that you have [conda](https://conda.io) or [mamba](https://mamba.readthedocs.io/) installed on your machine. We recommend `mamba`, as it's a faster drop-in replacement for `conda`. Using either one, you can create the environment:

        # using mamba
        mamba env create -f environment.yaml
        conda activate euro-calliope
        snakemake --use-conda --list # test your installation

        # using conda
        conda env create -f environment.yaml
        conda activate euro-calliope
        snakemake --use-conda --conda-frontend conda --list # test your installation

3. Install a Gurobi license on your computer ([academic license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) comes at no cost), or [choose a different solver](./customisation.md#manual).

4. Create an account at the Copernicus Climate Data Service and a `$HOME/.cdsapirc` file with your credentials; see their [How To](https://cds.climate.copernicus.eu/api-how-to) (you do not need to manually install the client).

## Build

By default, the entire model will be built in the `./build/model` folder when you execute the workflow. Execute the workflow like so:

```bash
snakemake --use-conda --cores <N_CORES>
```

`N_CORES` is the number of cores of your machine you want to use. It can be anything between `1` and `all`. Please have a look at [Snakemake's documentation](https://snakemake.readthedocs.io) for more information.

## Run

The build step creates all individual components of `euro-calliope`, like technologies and time series. These can be combined to eventually build a final model to run simulations with. For an example of such a model, see `./build/model/{resolution}/example-model.yaml`. It is a complete Calliope model and can be used like any other, for example like so:

```bash
calliope run ./build/model/national/example-model.yaml
```

For more information on how to use Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/).

## Customise

Without any modifications, the workflow will rebuild the [pre-built models](./pre-built.md). These models are examples and very likely require customisation to fit your purpose. Once you've managed to build and run them, it's a good point in time to learn about [customisation options in euro-calliope](./customisation.md).

## The repository

The following is an overview over the most important folders in the repository and their content.

* `build/model`: Contains the entire model after the build step, including Calliope definition files and data (does not exist initially).
* `config`: Files with configuration parameters that customise the model build process.
* `data/automatic`: Contains data automatically downloaded from third-party sources within the model build process (does not exist initially).
* `data/euro-calliope-dataset`: Contains data automatically downloaded from the [Euro-Calliope datasets repository](https://github.com/calliope-project/euro-calliope-datasets) within the model build process (does not exist initially).
* `docs`: Documentation of the model and the build process.
* `envs`: Files defining the conda environments that are used to build the model.
* `lib`: Library code in form of a Python package that is reused in many places of this repository.
* `rules`: Snakemake files defining the build process.
* `scripts`: Scripts that preprocess data and build the model components.
* `templates`: Templates of Calliope model files.
* `tests`: Tests of the models, scripts, and library code.

## On a cluster

### Build on a cluster

You may want to build the model on a cluster. While you can build euro-calliope on [any cluster that is supported by Snakemake](https://snakemake.readthedocs.io/en/v6.1.1/executing/cluster.html), our default configuration is targeted at, and tested on, ETH's Euler cluster. To build the model on Euler, use the following command:

```bash
snakemake --use-conda --profile config/euler
```

If you want to run on another cluster, read [snakemake's documentation on cluster execution](https://snakemake.readthedocs.io/en/v6.1.1/executable.html#cluster-execution) and take `config/euler` as a starting point.

### Work local, build on remote

If you are like us, you may want to work locally (to change configuration parameters, add modules etc), but execute remotely on the cluster. We support this workflow through three Snakemake rules: `send`, `receive`, and `clean_cluster_results`. It works like the following.

First, start local and make sure the `cluster-sync` configuration parameters fit your environment. Next, run `snakemake --use-conda --cores 1 send` to send the entire repository to your cluster. On the cluster, execute the workflow with Snakemake ([see above](./build.md#build-on-a-cluster)). After the workflow has finished, download results by locally running `snakemake --use-conda --cores 1 receive`. By default, this will download results into `build/cluster`.

This workflow works iteratively too. After analysing your cluster results locally, you may want to make changes locally, send these changes to the cluster (`snakemake --use-conda --cores 1 send`), rerun on the cluster, and download updated results (`snakemake --use-conda --cores 1 receive`).

To remove cluster results on your local machine, run `snakemake --use-conda --cores 1 clean_cluster_results`.

### Be notified of build successes or fails

 As the execution of this workflow may take a while, you can be notified whenever the execution terminates either successfully or unsuccessfully. Notifications are sent by email. To activate notifications, add the email address of the recipient to the configuration key `email`. You can add the key to your configuration file, or you can run the workflow the following way to receive notifications:

```bash
snakemake --use-conda --config email=<your-email>
```
