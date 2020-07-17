"""Applies config parameters to template files."""
from pathlib import Path

import jinja2

import filters


def parameterise_template(path_to_template, path_to_biofuel_costs, scaling_factors,
                          capacity_factors, biofuel_efficiency, max_power_density,
                          path_to_result):
    """Applies config parameters to template files."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    with open(path_to_biofuel_costs, "r") as f_biofuel_costs:
        biofuel_fuel_cost = float(f_biofuel_costs.readline())

    path_to_template = Path(path_to_template)
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(path_to_template.parent))
    env.filters['unit'] = filters.unit
    rendered = env.get_template(path_to_template.name).render(
        scaling_factors=scaling_factors,
        capacity_factors=capacity_factors,
        max_power_density=max_power_density,
        biofuel_fuel_cost=biofuel_fuel_cost,
        biofuel_efficiency=biofuel_efficiency
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        scaling_factors=snakemake.params["scaling_factors"],
        capacity_factors=snakemake.params["capacity_factors"],
        max_power_density=snakemake.params["max_power_density"],
        biofuel_efficiency=snakemake.params["biofuel_efficiency"],
        path_to_biofuel_costs=snakemake.input.biofuel_cost,
        path_to_result=snakemake.output[0]
    )
