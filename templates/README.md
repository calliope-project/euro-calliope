# Pre-built Euro-Calliope

Ready to use models of the European electricity system built for use in _Calliope_. Models are available on three different spatial resolutions: continental, national, and regional.

In addition, Euro-Calliope models can be built manually which adds more configuration options. To learn how to build Euro-Calliope manually, head over to [Euro-Calliope's documentation](https://euro-calliope.readthedocs.io).

## At a glance

Euro-Calliope models the European electricity system with each location representing an administrative unit. It is built on three spatial resolutions: on the continental level as a single location, on the national level with 34 locations, and on the regional level with 497 locations. At each location, renewable generation capacities (wind, solar, bioenergy) and balancing capacities (battery, hydrogen) can be built. In addition, hydro electricity and pumped hydro storage capacities can be built up to the extent to which they exist today. All capacities are used to satisfy electricity demand on all locations where demand is based on historic data. Locations are connected through transmission lines of either unrestricted capacity or projections. Using [Calliope](https://www.callio.pe), the model is formulated as a linear optimisation problem with total monetary cost of all capacities as the minimisation objective. Due to the flexibility of Calliope and the availability of the routines building the model all components can be adapted to the modeller's needs.

## Prepare

1. Install a Gurobi license on your computer ([academic license](https://www.gurobi.com/downloads/end-user-license-agreement-academic/) comes at no cost), or [choose a different solver](https://euro-calliope.readthedocs.io/en/latest/model/customisation/#manual-changes).

2. Install Calliope and all required dependencies. The easiest way to do so is using [conda](https://conda.io/) or [mamba](https://mamba.readthedocs.io/). Using conda, you can install Calliope:

```bash
cd pre-built-euro-calliope-v1.2.0
conda env create -f environment.yaml
conda activate euro-calliope
```

## Run

There are three models in the directory of the pre-builts -- one for each of the three spatial resolutions continental, national, and regional. You can run all three models out-of-the-box, but you may want to modify the model. By default, the model runs for the first day of January only. To run the example model on the continental resolution type:

```bash
calliope run ./continental/example-model.yaml
```

## Customise

The pre-built models are examples and very likely require customisation to fit your purpose. Once you've managed to run them, it's a good point in time to learn about [model customisation options in Euro-Calliope](https://euro-calliope.readthedocs.io/en/latest/model/customisation/).

## More information

For more information on Euro-Calliope and how to use and modify the models, see [Euro-Calliope's documentation](https://euro-calliope.readthedocs.io).

## License and attribution

Euro-Calliope is developed and maintained within the [Calliope project](https://www.callio.pe).

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

Contains modified Copernicus Atmosphere Monitoring Service information 2020. Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains.

Contains modified data from [Renewables.ninja](https://www.renewables.ninja/).

Contains modified data from [Open Power System Data](https://open-power-system-data.org).



