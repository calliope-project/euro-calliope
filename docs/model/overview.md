# Model overview

No matter whether you have downloaded Euro-Calliope's pre-builts or you have built the models yourself successfully, you will have a set of files in front of you. Let's make sense of these files.

## File structure

By default, Euro-Calliope is a set of three models on different spatial resolutions: continental, national, and regional.
All files required to run each resolution-specific model are within subfolders named by the resolution.
All files in the root directory are independent of the model configuration or data; they are not necessary to run a model.
Within each resolution-specific model directory, there is a subdirectory for technology definitions (stored in human-readable YAML files) and a subdirectory for timeseries data (stored in CSV files)

```
├── {resolution}                                            <- An individual folder for each spatial resolution.
│   ├── timeseries                                          <- All timeseries data CSV files.
|   |   |── supply
|   |   |   └── capacityfactors-{technology}.csv            <- Timeseries of capacityfactors of all renewables.
|   |   └── demand
|   |   |   └── electricity-demand.csv                      <- Timeseries of electricity demand on each node.
|   ├── techs                                               <- All technology definition YAML files.
|   |   |── {technology-class}                              <- Calliope base technology classes (one of `supply`, `demand`, `storage`, `transmission`).
|   |   |   └── {technology-class}-{technology-group}.csv   <- Definition of a technology (or group of technologies) relevant to the base technology, and the allocation of that technology to nodes in the model.
│   ├── example-model.yaml                                  <- Calliope model definition.
|   ├── interest-rate.yaml                                  <- Interest rates of all capacity investments.
|   ├── scenarios.yaml                                      <- scenario names which can be used to override the base model configuration.
│   └── locations.yaml                                      <- Defines all nodes in the model, including the coordinates defining their centroids.
├── build-metadata.yaml                                     <- Metadata of the build process.
├── environment.yaml                                        <- Conda file defining an environment to run the model in.
└── README.md                                               <- Basic documentation (pre-builts only).
```

## Units of quantities

The units of quantities within the models are the following:

* power: 100 GW
* energy: 100 GWh
* area: 10,000 km²
* monetary cost: billion EUR

All data going into Calliope and all Calliope result data will be given using these units.
While they may be unusual, these units lead to a numerical model that is well suited for the interior-point solution algorithm that is used by default.
The units are tuned so as to work best for models with a time resolution of a few hours and a duration of one year.
For other types of problems, or other solution algorithms, the units may need to be changed to avoid numerical issues within the solver.

When you run the workflow, you can easily change the units and scale all values using the `scaling-factor` configuration parameters, see [Configuration](./customisation.md#configuration).
The base units on which the scaling factors are applied are `1 MW`, `1 MWh`, `EUR`, and `km2`.
So for example, the default unit for energy (100 GWh) is derived by scaling the base unit (1 MWh) with a scaling factor of `0.00001`.

## Components and assumptions

For an in-depth description of all model components and the data preprocessing steps, please read the [open-access article introducing Euro-Calliope in Joule](https://doi.org/10.1016/j.joule.2020.07.018).

## Your feedback

Do you have an idea how we could improve Euro-Calliope? Did you find an issue with the model or its execution? [Please let us know by contributing to our issue tracker.](https://github.com/calliope-project/euro-calliope/issues/new/choose)
