# Customising your model

The example models that come with Euro-Calliope are complete models that you can run using Calliope.
However, it's unlikely that these models fit your purpose optimally.
To make them fit-for-purpose, Euro-Calliope models can be configured, adapted, and extended.
You have the following three options:

1. [Manual changes](./customisation.md#manual-changes)
2. [Importing modules](./customisation.md#importing-modules)
3. [Overrides and scenarios](./customisation.md#overrides-and-scenarios)
4. [Rebuild model](./customisation.md#rebuild)

## Manual changes

With the Calliope model in your hands, you will be able to change any model parameter, any technology specifics, and the model definition to your liking.
This kind of customisation can be useful to get to know the model and its parameters.
To create reliable results, we do advise making manual changes only to the model definition (`example-model.yaml`) as this makes it possible to trace those changes later.
A typical customisation here would be to change the solver from `gurobi` to an open-source solver, e.g. `cbc` (see [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/user/config_defaults.html#run-configuration)).
We consider all Euro-Calliope model subcomponents (everything other than the model definition itself) as a toolbox from which you can choose to define your model -- see the [Import customisation option](./customisation.md#importing-modules).

## Importing modules

The `example-model.yaml` definition file in each resolution sub-directory (e.g. `national/example-model.yaml`) specifies a list of other files to bring together to describe the model (under the `import` key).
This list can be changed by the modeller to select a combination of different files (see also [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#files-that-define-a-model)).
These files represent "modules" of the model definition and contain everything necessary for a given technology or technology group to exist.
For instance, `techs/supply/hydro.yaml` defines two technologies (under the `techs` key) which will convert river flows into electricity.
It also places that technology in every relevant modelled location (under the `locations` key), along with any location-specific information that is needed; in this case, the maximum capacity of hydropower in that location.
Finally, there are potential overrides defined, which is an additional layer of customisation described further in the [overrides section](./customisation.md#overrides-and-scenarios).
Only by including `techs/supply/hydro.yaml` in the list of imports in `example-model.yaml` will the defined hydropower technologies exist in the built Calliope model.
By default, the example model definition imports all modules except electricity transmission, so you can simply remove any modules from the list of imports if you do not want to consider that technology / technology group in your study.

The modules are in the `techs` subdirectory of each spatial resolution (e.g. `national/techs/...`).
Here, we describe each module in terms of the technologies they contain (`calliope name`:`full name`) and the overrides they make available (`calliope name`: `override description`) .

??? note "demand/electricity.yaml"

    === "Technologies"

        **demand_elec**: Electricity demand

??? note "storage/electricity.yaml"

    === "Technologies"

        **battery**: Battery storage

        **hydrogen**: Hydrogen power storage

    === "Overrides"

        **exclusive-energy-to-power-ratios**: Constrain the energy to power ratios of battery and hydrogen storage in a way that they do not overlap (in Calliope terms, energy="storage capacity", power="energy capacity").
        Battery storage is constrained to a ratio of ≤4h while hydrogen is constrained to a ratio of ≥4h.
        The ratio is derived from typical values of commercial lithium-ion batteries available today (2021).
        Constraining hydrogen storage as well ensures it does not directly compete with battery storage, but is used instead for durations of fours hours and longer.


??? note "storage/hydro.yaml"

    === "Technologies"

        **pumped_hydro**: Pumped hydro power storage

    === "Overrides"

        **no-hydro-storage-fixed-cost**: Set installation costs of pumped hydro storage to zero.

        **freeze-hydro-storage-capacities**: Force discharge and storage capacities of hydropower to today's levels instead of using today's level as the maximum (you may also want to apply "no-hydro-storage-fixed-cost").

??? note "supply/biofuel.yaml"

    === "Technologies"

        **biofuel**: Biofuel

??? note "supply/hydro.yaml"

    === "Technologies"

        **hydro_reservoir**: Hydro electricity with a reservoir.

        **hydro_run_of_river**: Run of river hydro electricity

    === "Overrides"

        **no-hydro-supply-fixed-cost**: Set installation costs of hydropower supply technologies to zero.

        **schroeder-hydro-cost**: Override hydropower supply cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014 with those from Schröder et al. (2013).

        **freeze-hydro-supply-capacities**: Force discharge and storage capacities of hydropower to today's levels instead of using today's level as the maximum (you may also want to apply "no-hydro-supply-fixed-cost").

??? note "supply/load-shedding.yaml"

    === "Technologies"

        **load_shedding**: Load shedding as last resort

    === "Overrides"

        **load-shedding**: Add an option to shed load at each location.
        You can use this to model blackouts, brownouts, or controlled shedding of load as a form of demand response.
        In Euro-Calliope, we model load shedding not as actual reduction of demand but as an unconstrained supply of electricity.
        This supply has high variable cost (see the load-shedding.yaml module) and no fixed cost.
        Due to its high cost, it will only be used when no other, less costly, option is available.


??? note "supply/nuclear.yaml"

    === "Technologies"

        **nuclear**: Nuclear power

??? note "supply/open-field-solar-and-wind-onshore.yaml"

    === "Technologies"

        **open_field_pv**: Open field PV

        **wind_onshore_competing**: Onshore wind competing with open field PV on land

        **wind_onshore_monopoly**: Onshore wind without land competition

    === "Overrides"

        **dea-renewable-cost-pv-open-field**: Override open field PV cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014 with those from the Danish Energy Agency.

        **dea-renewable-cost-wind-onshore**: Override onshore wind cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014 with those from the Danish Energy Agency.

??? note "supply/rooftop-solar.yaml"

    === "Technologies"

        **roof_mounted_pv**: Roof mounted PV

    === "Overrides"

        **dea-renewable-cost-pv-roof-mounted**: Override cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014 with those from the Danish Energy Agency

        **directional-rooftop-pv**: By default, Euro-Calliope contains a single technology for rooftop PV.
        This technology comprises the total rooftop PV potential in each location, in particular including east-, west-, and north-facing rooftops.
        While this allows you the model fully exploit the potential of rooftop PV, it leads to less than optimal capacity factors as long as the potential is not fully exploited.
        That is because in reality, one would likely first exploit all south-facing rooftops, then east- and west-facing rooftops, and only then -- if at all -- north-facing rooftops.
        When using this override, there are three technologies instead of just one for rooftop PV.
        The three technologies comprise (1) south-facing PV (on either south-facing or flat rooftops), (2) east- and west-facing PV, and (3) north-facing PV.
        This leads to higher capacity factors of rooftop PV as long as the potential of rooftop PV is not fully exploited.
        However, this also increases the complexity of the model.

??? note "supply/wind-offshore.yaml"

    === "Technologies"

        **wind_offshore**: Offshore wind

    === "Overrides"

        **dea-renewable-cost-wind-offshore**: Override offshore wind cost and lifetime projections from the JRC Energy Technology Reference Indicator 2014 with those from the Danish Energy Agency.

??? note "transmission/electricity-entsoe.yaml"

    === "Technologies"

        **ac_transmission**: High voltage AC transmission line

        **free_transmission**: Local power transmission

??? note "transmission/electricity-linked-neighbours.yaml"

    === "Technologies"

        **ac_transmission**: High voltage AC transmission line

        **free_transmission**: Local power transmission


## Overrides and scenarios

Calliope [overrides](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#scenarios-and-overrides) enable models to be easily manipulated.
An override named `freeze-hydro-supply-capacities` can be used for example in this way:

``` bash
calliope run build/models/continental/example-model.yaml --scenario=freeze-hydro-supply-capacities
```

``` python
import calliope
model = calliope.Model("build/models/continental/example-model.yaml", scenario="freeze-hydro-supply-capacities")
model.run()
```

Overrides can also be chained, enabling multiple scenarios to built up from multiple overrides.
For instance, `freeze-hydro-supply-capacities` and `freeze-hydro-storage-capacities` can be combined to `freeze-hydro-supply-capacities,freeze-hydro-storage-capacities`.

You can also define your own overrides to manipulate any model component.
We recommend you add these overrides into the model definition YAML file, to ensure they are easy to trace.

In Calliope, [scenarios](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#scenarios-and-overrides) are groups of overrides and/or other scenarios.
In Euro-Calliope, it can be helpful to define scenarios to help group similar overrides together.
For instance, cost overrides from the Danish Energy Agency are defined in various files, since they are loaded in alongside the technologies they affect (the option to override offshore wind costs only exists when you load the `techs/supply/wind-offshore.yaml` module).
You can pre-define scenarios in your model definition file, such as:

``` yaml
scenarios:
    dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-wind-offshore, dea-renewable-cost-pv-roof-mounted]
```

Then you can load in the scenario into calliope as follows:

``` bash
calliope run build/models/continental/example-model.yaml --scenario=dea-renewable-cost
```

In the above example, if you choose not to load the offshore wind module into your model, then the scenario you define would become:

``` yaml
scenarios:
    dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-pv-roof-mounted]
```

Similar overrides which you may wish to group together are:

``` yaml
freeze-hydro-capacities: [freeze-hydro-supply-capacities, freeze-hydro-storage-capacities]
dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-wind-offshore, dea-renewable-cost-pv-roof-mounted]
no-hydro-fixed-cost: [no-hydro-supply-fixed-cost, no-hydro-storage-fixed-cost]
```

## Rebuild

When above options do not provide enough flexibility for you, you can rebuild the model using Euro-Calliope's workflow and use the [customisation options of the workflow](../workflow/customisation.md).
