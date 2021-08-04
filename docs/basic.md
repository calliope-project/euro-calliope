# Basic usage and the pre-built models

The pre-built models of euro-calliope are regular Calliope models that you can use like any other Calliope model. If you are unfamiliar with Calliope, we'll show you below how to run the model. After going through these first steps, we advice you to head over to [Calliope's documentation](https://calliope.readthedocs.io/en/stable/) to learn its basic usage.

## Prepare

1. [Download](https://doi.org/10.5281/zenodo.3949552) the pre-built models.

2. You need a Gurobi license installed on your computer. You may as well choose another solver than Gurobi. See [Calliope's documentation](https://calliope.readthedocs.io/en/stable/user/config_defaults.html?highlight=solver#run-configuration) to understand how to switch to another solver.

2. Install Calliope and all required dependencies. The easiest way to do so is using [conda](https://conda.io/) or [mamba](https://mamba.readthedocs.io/). Using conda, you can install Calliope:

```bash
cd pre-built-euro-calliope-v1.1
conda env create -f environment.yaml
conda activate euro-calliope
```

## Run

There are three models in the directory of the pre-builts -- one for each of the three spatial resolutions continental, national, and regional. You can run all three models out-of-the-box, but you may want to modify the model. By default, the model runs for the first day of January only. To run the example model on the continental resolution type:

```bash
calliope run ./continental/example-model.yaml
```

For more information on how to use and modify Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io) and [euro-calliope's customisation options](./customisation.md).
