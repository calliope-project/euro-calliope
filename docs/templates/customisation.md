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
We consider all Euro-Calliope model subcomponents (everything other than the model definition itself) as a toolbox from which you can choose to define your model -- see the [Import customisation option](./customisation.md#imports).

## Importing modules

The `example-model.yaml` definition file in each resolution sub-directory (e.g. `national/example-model.yaml`) specifies a list of other files to bring together to describe the model (under the `import` key).
This list can be changed by the modeller to select a combination of different files (see also [Calliope's documentation](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#files-that-define-a-model)).
These files represent "modules" of the model definition and contain everything necessary for a given technology or technology group to exist.
For instance, `techs/supply/hydro.yaml` defines two technologies (under the `techs` key) which will convert river flows into electricity.
It also places that technology in every relevant modelled location (under the `locations` key), along with any location-specific information that is needed; in this case, the maximum capacity of hydropower in that location.
Finally, there are potential overrides defined, which is an additional layer of customisation described further in the [overrides section](./customisation.md#overrides).
Only by including `techs/supply/hydro.yaml` in the list of imports in `example-model.yaml` will the defined hydropower technologies exist in the built Calliope model.
By default, the example model definition imports all modules except electricity transmission, so you can simply remove any modules from the list of imports if you do not want to consider that technology / technology group in your study.

The modules are in the `techs` subdirectory of each spatial resolution (e.g. `national/techs/...`).
Here, we describe each module in terms of the technologies they contain (`calliope name`:`full name`) and the overrides they make available (`calliope name`: `override description`) .

{% for module_name, module_contents in modules.items() %}
<details>
<summary><code>{{ module_name }}</code></summary>
<details>
<summary><i>Technologies</i></summary>
<ul>
{% for tech, tech_name in module_contents.techs.items() %}
<li><b>{{ tech }}</b>: {{ tech_name }}</li>
{% endfor %}
</ul>
</details>
{% if module_contents.overrides is defined %}
<details>
<summary><i>Overrides</i></summary>
<ul>
{% for override_name, override in module_contents.overrides.items() %}
<li><b>{{ override_name }}</b>: {{ override }}</li>
{% endfor %}
</ul>
</details>
{% endif %}
</details>
{% endfor %}
<br>

## Overrides and scenarios

Calliope [overrides](https://calliope.readthedocs.io/en/v0.6.7/user/building.html#scenarios-and-overrides) enable models to be easily manipulated.
An override named `freeze-hydro-supply-capacities` can be used for example in this way:

```bash
calliope run build/models/continental/example-model.yaml --scenario=freeze-hydro-supply-capacities
```

```python
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

```yaml
scenarios:
    dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-wind-offshore, dea-renewable-cost-pv-roof-mounted]
```

Then you can load in the scenario into calliope as follows:

```bash
calliope run build/models/continental/example-model.yaml --scenario=dea-renewable-cost
```

In the above example, if you choose not to load the offshore wind module into your model, then the scenario you define would become:

```yaml
scenarios:
    dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-pv-roof-mounted]
```

Similar overrides which you may wish to group together are:

```yaml
freeze-hydro-capacities: [freeze-hydro-supply-capacities, freeze-hydro-storage-capacities]
dea-renewable-cost: [dea-renewable-cost-pv-open-field, dea-renewable-cost-wind-onshore, dea-renewable-cost-wind-offshore, dea-renewable-cost-pv-roof-mounted]
no-hydro-fixed-cost: [no-hydro-supply-fixed-cost, no-hydro-storage-fixed-cost]
```

## Rebuild

When above options do not provide enough flexibility for you, you can rebuild the model using Euro-Calliope's workflow and use the [customisation options of the workflow](../workflow/customisation.md).
