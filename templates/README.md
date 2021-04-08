# pre-built euro-calliope

Ready to use models of the European electricity system built using Calliope. Models are available on three different spatial resolutions: continental, national, and regional.

In addition, euro-calliope models can be built manually which adds more configuration options. To build euro-calliope manually, head over to [GitHub](https://github.com/calliope-project/euro-calliope).

## At a glance

euro-calliope models the European electricity system with each location representing an administrative unit. It is built on three spatial resolutions: on the continental level as a single location, on the national level with 34 locations, and on the regional level with 497 locations. On each node, renewable generation capacities (wind, solar, bioenergy) and balancing capacities (battery, hydrogen) can be built. In addition, hydro electricity and pumped hydro storage capacities can be built up to the extent to which they exist today. All capacities are used to satisfy electricity demand on all locations which is based on historic data. Locations are connected through transmission lines of unrestricted capacity. Using [Calliope](https://www.callio.pe), the model is formulated as a linear optimisation problem with total monetary cost of all capacities as the minimisation objective. The pre-built models can be manipulated by updating any of the files. In addition to the pre-built models, models can be built manually. Manual builds provide more flexibility in adapting and configuring the model. To build euro-calliope manually, head over to [GitHub](https://github.com/calliope-project/euro-calliope).

## Get ready to run the models

1. You need a Gurobi license installed on your computer. You may as well choose another solver than Gurobi. See [Calliope's documentation](https://calliope.readthedocs.io/en/stable/user/config_defaults.html?highlight=solver#run-configuration) to understand how to switch to another solver.

2. You need to have Calliope and Gurobi installed in your environment. The easiest way to do so is using [conda](https://conda.io/docs/index.html). Using conda, you can create a conda environment from within you can build the model:

```
conda env create -f environment.yaml
conda activate euro-calliope
```

## Run the models

There are three models in this directory -- one for each of the three spatial resolutions continental, national, and regional. You can run all three models out-of-the-box, but you may want to modify the model. By default, the model runs for the first day of January only. To run the example model on the continental resolution type:

```Bash
$ calliope run ./continental/example-model.yaml
```

For more information on how to use and modify Calliope models, see [Calliope's documentation](https://calliope.readthedocs.io).

## Manipulating the model using overrides

Calliope [overrides](https://calliope.readthedocs.io/en/stable/user/building.html#scenarios-and-overrides) allow to easily manipulate models. An override named `freeze-hydro-capacities` can be used for example in this way:

```bash
calliope run build/model/continental/example-model.yaml --scenario=freeze-hydro-capacities
```

You can define your own overrides to manipulate any model component. The following overrides are built into euro-calliope:

> directional-rooftop-pv

By default, euro-calliope contains a single technology for rooftop PV. This technology comprises the total rooftop PV potential in each location, in particular including east-, west-, and north-facing rooftops. While this allows to fully exploit the potential of rooftop PV, it leads to less than optimal capacity factors as long as the potential is not fully exploited. That is because, one would likely first exploit all south-facing rooftop, then east- and west-facing rooftops, and only then -- if at all -- north-facing rooftops. By default, euro-calliope cannot model that.

When using the `directional-rooftop-pv` override, there are three instead of just one technologies for rooftop PV. The three technologies comprise (1) south-facing and flat rooftops, (2) east- and west-facing rooftops, and (3) north-facing rooftops. This leads to higher capacity factors of rooftop PV as long as the potential of rooftop PV is not fully exploited. However, this also increases the complexity of the model.

> freeze-hydro-capacities

By default, euro-calliope allows capacities of run-of-river hydro, reservoir hydro, and pumped storage hydro capacities up to today's levels. Alternatively, it's possible to freeze these capacities to today's levels using the `freeze-hydro-capacities` override.

## Model components

The models contain the following files. All files in the root directory are independent of the spatial resolution. All files that depend on the spatial resolution are within subfolders named by the resolution.

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
├── environment.yaml                       <- Conda file defining an environment to run the model in.
├── interest-rate.yaml                     <- Interest rates of all capacities.
├── link-techs.yaml                        <- Definition of link technologies.
├── README.md                              <- The file you are currently looking at.
├── renewable-techs.yaml                   <- Definition of supply technologies.
└── storage-techs.yaml                     <- Definition of storage technologies.
```

## Units of quantities

The units of quantities within the models are the following:

* power: {{ (1 / scaling_factors.power) | unit("MW", parenthesis=False) }}
* energy: {{ (1 / scaling_factors.power) | unit("MWh", parenthesis=False) }}
* area: {{ (1 / scaling_factors.area) | unit("km2", parenthesis=False) }}
* monetary cost: {{ (1 / scaling_factors.monetary) | unit("EUR", parenthesis=False) }}

These units were chosen in order to minimise numerical issues within the optimisation algorithm.

## License and attribution

euro-calliope has been developed and is maintained by Tim Tröndle, IASS Potsdam.

<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

Contains modified Copernicus Atmosphere Monitoring Service information 2020. Neither the European Commission nor ECMWF is responsible for any use that may be made of the Copernicus information or data it contains.

Contains modified data from [Renewables.ninja](https://www.renewables.ninja/).

Contains modified data from [Open Power System Data](https://open-power-system-data.org).



