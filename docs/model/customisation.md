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
To create reliable results, we advise making manual changes only to the model definition (`example-model.yaml`) as this makes it possible to trace those changes later.
A typical customisation here would be to change the solver from `gurobi` to an open-source solver, e.g. `cbc` (see [Calliope's documentation](https://calliope.readthedocs.io/en/v{{ calliope_version }}/user/config_defaults.html#run-configuration)).
We consider all Euro-Calliope model subcomponents (everything other than the model definition itself) as a toolbox from which you can choose to define your model -- see the [Import customisation option](./customisation.md#importing-modules).

## Importing modules

The `example-model.yaml` definition file in each resolution sub-directory (e.g. `national/example-model.yaml`) specifies a list of other files to bring together to describe the model (under the `import` key).
This list can be changed by the modeller to select a combination of different files (see also [Calliope's documentation](https://calliope.readthedocs.io/en/v{{ calliope_version }}/user/building.html#files-that-define-a-model)).
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

??? note "demand/electrified-transport.yaml"

    === "Technologies"

        **demand_road_transport_electrified_uncontrolled**: Share of electrified road transport demand which is uncontrolled.

        **demand_road_transport_historic_electrified_uncontrolled**: Removes historically electrified road transport demand to avoid double counting. It is assumed uncontrolled.

        **demand_road_transport_electrified_controlled**: Share of electrified road transport demand whose charging is optimised by the solver.

    === "Overrides"

        **keep-historic-electricity-demand-from-road-transport**: Keep historically electrified road transport demand. Historically electrified road transport demand is deleted by default, as it is already considered in historic electricity demand and would thus be counted twice. Using this override together with Euro-Calliope's default electricity demand is not advised.

        **(year)_transport_controlled_electrified_demand**: Total electrified road transport demand whose charging is optimised by the solver.

??? note "demand/heat.yaml"

    === "Technologies"

        **demand_heat**: Combined space heat and hot water demand.

??? note "storage/electricity.yaml"

    === "Technologies"

        **battery**: Battery storage

        **hydrogen**: Hydrogen power storage

    === "Overrides"

        **exclusive-energy-to-power-ratios**: Constrain the energy to power ratios of battery and hydrogen storage in a way that they do not overlap (in Calliope terms, energy="storage capacity", power="energy capacity").
        Battery storage is constrained to a ratio of ≤4h while hydrogen is constrained to a ratio of ≥4h.
        The ratio is derived from typical values of commercial lithium-ion batteries available today (2021).
        Constraining hydrogen storage as well ensures it does not directly compete with battery storage, but is used instead for durations of fours hours and longer.

??? note "storage/heat.yaml"

    === "Technologies"

        **heat_storage_small**: Abstract [technology group](https://calliope.readthedocs.io/en/v0.6.10/user/advanced_features.html#using-tech-groups-to-group-configuration).
        This "technology" only becomes part of the model when defining technologies in the overrides of this file.

        **hp_heat_storage_small**: Storage buffer for heat pumps which inherits from the `heat_storage_small` abstract technology group, assuming a domestic (small scale) application.

        **electric_heater_heat_storage_small**: Storage buffer for direct electric heaters which inherits from the `heat_storage_small` abstract technology group, assuming a domestic (small scale) application.

        **biofuel_heat_storage_small**: Storage buffer for biofuel boilers which inherits from the `heat_storage_small` abstract technology group, assuming a domestic (small scale) application.

        **methane_heat_storage_small**: Storage buffer for methane boilers which inherits from the `heat_storage_small` abstract technology group, assuming a domestic (small scale) application.

??? note "storage/hydro.yaml"

    === "Technologies"

        **pumped_hydro**: Pumped hydro power storage

    === "Overrides"

        **no-hydro-storage-fixed-cost**: Set installation costs of pumped hydro storage to zero.

        **freeze-hydro-storage-capacities**: Force discharge and storage capacities of hydropower to today's levels instead of using today's level as the maximum (you may also want to apply "no-hydro-storage-fixed-cost").

??? note "supply/biofuel.yaml"

    === "Technologies"

        **biofuel**: Biofuel supply, limited per model location to an hourly flow that can be stored before use in downstream technologies.

??? note "supply/electrified-biofuel.yaml"

    === "Technologies"

        **electrified_biofuel**: Electrified biofuel supply (assuming anaerobic digestion).
        This should be used in an electricity-only model, where biofuel supply is assumed to only be used for direct generation of electricity.
        This simplifies the model compared to using `supply/biofuel.yaml` and `conversion/electricity-from-biofuel.yaml` by not introducing the `biofuel` energy carrier.

??? note "conversion/heat-from-biofuel.yaml"

    === "Technologies"

        **biofuel_boiler**: Biofuel-consuming boiler.

        **biofuel_tech_heat_to_demand**: "Dummy" technology to convert biofuel boiler output to a carrier that can be used to meet heat demand.

??? note "conversion/electricity-from-biofuel.yaml"

    === "Technologies"

        **electricity_from_biofuel**: Biofuel-consuming electricity production facility (assuming anaerobic digestion).

??? note "conversion/heat-from-electricity.yaml"

    === "Technologies"

        **heat_pump**: Heat pump.

        **heat_pump_tech_heat_to_demand**: Dummy technology to convert heat pump output to a carrier that can be used to meet heat demand.

        **electric_heater**: Direct electric heater.

        **electric_heater_tech_heat_to_demand**: Dummy technology to convert electric heater output to a carrier that can be used to meet heat demand.

??? note "conversion/heat-from-methane.yaml"

    === "Technologies"

        **methane_boiler**: Natural gas / methane boiler

        **methane_tech_heat_to_demand**: "Dummy" technology to convert methane boiler output to a carrier that can be used to meet heat demand.

??? note "conversion/synfuels-from-hydrogen.yaml"

    === "Technologies"

        **hydrogen_to_liquids**: Hydrogen+CO<sub>2</sub> to liquid fuels (diesel & kerosene) converter.

        **hydrogen_to_methanol**: Hydrogen+CO<sub>2</sub> to methanol converter.

        **hydrogen_to_methane**: Hydrogen+CO<sub>2</sub> to methane converter.

        **dac**: Direct air CO<sub>2</sub> capture.

??? note "conversion/fuels-from-synfuels.yaml"

    === "Technologies"

        **syn_diesel_converter**: "Dummy" technology to convert _synthetic_ diesel to a carrier that can be used to meet diesel demand.

        **syn_methane_converter**: "Dummy" technology to convert _synthetic_ methane to a carrier that can be used to meet methane demand.

        **syn_kerosene_converter**: "Dummy" technology to convert _synthetic_ kerosene to a carrier that can be used to meet kerosene demand.

        **syn_methanol_converter**: "Dummy" technology to convert _synthetic_ methanol to a carrier that can be used to meet methanol demand.

??? note "conversion/synfuels-from-biofuel.yaml"

    === "Technologies"

        **biofuel_to_liquids**: Biofuel to liquid fuels (diesel & kerosene) converter.

        **biofuel_to_diesel**: Biofuel to vehicle fuel (assumed diesel) converter.

        **biofuel_to_methanol**: Biofuel to methanol converter.

        **biofuel_to_methane**: Biofuel to methane converter.

??? note "conversion/hydrogen-from-electricity.yaml"

    === "Technologies"

        **electrolysis**: Hydrogen by electrolysis, assuming an equal combination from the primary types of electrolysers: SOEC, PEM, and Alkaline.

??? note "supply/historic-electrified-heat.yaml"

    === "Technologies"

        **demand_heat_historic_electrified**: Removes historically electrified heat demand to avoid double counting in the electricity demand profile.

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
        In Euro-Calliope, we model load shedding not as actual reduction of demand but as an unconstrained supply of electricity.
        This supply has high variable cost (see `tech-cost.yaml` parameter file) and no fixed cost.
        Due to its high cost, it will only be used when no other, less costly, option is available.

        Calliope provides a built-in mechanism that is similar: [`ensure-feasibility`](https://calliope.readthedocs.io/en/v{{ calliope_version }}/user/building.html#allowing-for-unmet-demand).
        The benefit of using the `load-shedding` override over Calliope's built-in mechanism is that it is more targeted towards modelling shedding of electrical load and provides more flexibility -- for example in terms of the cost of shed load.

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

??? note "transmission/synfuel.yaml"

    This file provides the ability to move synthetic fuels between any two locations without have explicit transmission links.
    To work as expected, an additional, custom constraint is required to equate annual export and import for each carrier at each location.

    === "Technologies"

        **syn_diesel_distribution_export**: Synthetic diesel exporter

        **syn_methane_distribution_export**: Synthetic methane exporter

        **syn_kerosene_distribution_export**: Synthetic kerosene exporter

        **syn_methanol_distribution_export**: Synthetic methanol exporter

        **syn_diesel_distribution_import**: Synthetic diesel importer

        **syn_methane_distribution_import**: Synthetic methane importer

        **syn_kerosene_distribution_import**: Synthetic kerosene importer

        **syn_methanol_distribution_import**: Synthetic methanol importer


## Overrides and scenarios

Calliope [overrides](https://calliope.readthedocs.io/en/v{{ calliope_version }}/user/building.html#scenarios-and-overrides) enable models to be easily manipulated.
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

In Calliope, [scenarios](https://calliope.readthedocs.io/en/v{{ calliope_version }}/user/building.html#scenarios-and-overrides) are groups of overrides and/or other scenarios.
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
