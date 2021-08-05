# The pre-built models

If you want an easy entry into using and exploring euro-calliope, our pre-built models provide the jump start you need.
The pre-built models are ordinary Calliope models that you can use like any other.
If you are unfamiliar with Calliope, we'll show you below how to run the model.
After going through these first steps, we advise you to head over to [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/) to learn its basic usage.

## Prepare

1. [Download](https://doi.org/10.5281/zenodo.3949552) the pre-built models.

2. Install a Gurobi license on your computer ([academic license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) comes at no cost), or [choose a different solver](./customisation.md#manual).

3. Install Calliope and all required dependencies.
The easiest way to do so is using [conda](https://conda.io/) or [mamba](https://mamba.readthedocs.io/).
Using conda, you can install Calliope:

```bash
cd pre-built-euro-calliope-v1.1
conda env create -f environment.yaml
conda activate euro-calliope
```

## Run

There are three models in the directory of the pre-builts -- one for each of the three spatial resolutions continental, national, and regional.
You can run all three models out-of-the-box, but you may want to modify the model.
By default, the model runs for the first day of January only.
To run the example model on the continental resolution type:

```bash
calliope run ./continental/example-model.yaml
```

For more information on how to use and modify Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io) and [euro-calliope's customisation options](./customisation.md).

## Customise

The pre-built models are examples and very likely require customisation to fit your purpose.
Once you've managed to run them, it's a good point in time to learn about [customisation options in euro-calliope](./customisation.md).
