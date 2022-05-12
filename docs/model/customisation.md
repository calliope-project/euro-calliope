# Customising your model

The example models that come with Euro-Calliope are complete models that you can run using Calliope.
However, it's unlikely that these models fit your purpose optimally.
To make them fit-for-purpose, Euro-Calliope models can be configured, adapted, and extended.
You have the following three options:

1. [Manual changes](./customisation.md#manual-changes)
2. [Imports](./customisation.md#file-imports)
3. [Overrides](./customisation.md#overrides)
4. [Rebuild model](./customisation.md#rebuild)

## Manual changes

With the Calliope model in your hands, you will be able to change any model parameter, any technology specifics, and the model definition to your liking.
This kind of customisation can be useful to get to know the model and its parameters.
To create reliable results, we do not advice to make any manual changes to anything but the model definition as this may impact traceability of your results.
For the model definition (as in `./{resolution}/example-model.yaml`) we encourage you to do manual changes.
A typical customisation here would be to change the solver from `gurobi` to an open-source solver, e.g. `cbc` (see [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/user/config_defaults.html#run-configuration)).
We consider all Euro-Calliope model components but the model definition as a toolbox from which you can choose to define your model -- see the [Import customisation option](./customisation.md#imports).

## Imports

The `example-model.yaml` definition file in each resolution sub-directory (e.g. `national/example-model.yaml`) specifies a list of other files to bring together to describe the model.
This list can be changed by the modeller to select a combination of different files (see also [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#files-that-define-a-model)).

### Transmission links

Transmission links between locations in your model depend on the imported model file.
When you do not import any file, all locations are isolated (modelling full autarky).

For the national and regional resolutions with more than one location, you can import `./{resolution}/link-all-neighbours.yaml` which includes links between all neighbouring regions and a selection of pre-defined sub-sea links, but has no capacity limits.
This file is imported in the example models and the pre-builts.

At the national resolution, transmission links can be set based on an ENTSO-E ten-year development plan 2020 scenario (`national/entsoe-tyndp-links.yaml`).
The ENTSO-E links define all existing and planned international connections, including their predicted net transfer capacities (NTCs).

## Overrides

Calliope [overrides](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#scenarios-and-overrides) enable models to be easily manipulated.
An override named `dea-renewable-cost` can be used for example in this way:

```bash
calliope run build/models/continental/example-model.yaml --scenario=dea-renewable-cost
```

You can define your own overrides to manipulate any model component.
The following overrides are built into Euro-Calliope.

### Cost

By default, Euro-Calliope uses cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014.
The `dea-renewable-cost` override allows to use the projections from the Danish Energy Agency instead for solar PV and wind power and `schroeder-hydro-cost` provides another source for the hydropower assumptions.
Using the override `no-hydro-fixed-cost` allows to only consider variable and O&M costs for hydropower.
This may make sense especially in combination with the `freeze-hydro-capacities` override (see below).

### directional-rooftop-pv

By default, Euro-Calliope contains a single technology for rooftop PV.
This technology comprises the total rooftop PV potential in each location, in particular including east-, west-, and north-facing rooftops.
While this allows to fully exploit the potential of rooftop PV, it leads to less than optimal capacity factors as long as the potential is not fully exploited.
That is because in reality, one would likely first exploit all south-facing rooftops, then east- and west-facing rooftops, and only then -- if at all -- north-facing rooftops. By default, Euro-Calliope cannot model that.

When using the `directional-rooftop-pv` override, there are three instead of just one technologies for rooftop PV.
The three technologies comprise (1) south-facing PV (on either south-facing or flat rooftops), (2) east- and west-facing PV, and (3) north-facing PV.
This leads to higher capacity factors of rooftop PV as long as the potential of rooftop PV is not fully exploited.
However, this also increases the complexity of the model.

### exclusive-energy-to-power-ratios

Adding the `exclusive-energy-to-power-ratios` override constrains the energy to power ratios of battery and hydrogen storage in a way that they do not overlap (in Calliope terms: energy="storage capacity", power="energy capacity").
Battery storage is constrained to a ratio of ≤4h while hydrogen is constrained to a ratio of ≥4h.
The ratio is derived from typical values of commercial lithium-ion batteries available today (2021).
Constraining hydrogen storage as well ensures it does not directly compete with battery storage, but is used instead for durations of fours hours and longer.

### freeze-hydro-capacities

By default, Euro-Calliope allows capacities of run-of-river hydro, reservoir hydro, and pumped storage hydro capacities up to today's levels.
Alternatively, it's possible to freeze these capacities to today's levels using the `freeze-hydro-capacities` override.

### load-shedding

The `load-shedding` override adds an option to shed load at each location.
You can use this to model blackouts, brownouts, or controlled shedding of load as a form of demand response.

In Euro-Calliope, we model load shedding not as actual reduction of demand but as an unconstrained supply of electricity.
This supply has high variable cost (see `tech-cost.yaml` parameter file) and no fixed cost.
Due to its high cost, it will only be used when no other, less costly, option is available.

Calliope provides a built-in mechanism that is similar: [`ensure-feasibility`](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#allowing-for-unmet-demand).
The benefit of using the `load-shedding` override over Calliope's built-in mechanism is that it is more targeted towards modelling shedding of electrical load and provides more flexibility -- for example in terms of the cost of shed load.

## Rebuild

When above options do not provide enough flexibility for you, you can rebuild the model using Euro-Calliope's workflow and use the [customisation options of the workflow](../workflow/customisation.md).
