"""Applies config parameters to template files."""
import jinja2


def parameterise_template(path_to_template, path_to_biofuel_costs, scaling_factors,
                          biofuel_efficiency, max_power_density, path_to_result):
    """Applies config parameters to template files."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]
    with open(path_to_biofuel_costs, "r") as f_biofuel_costs:
        biofuel_fuel_cost = float(f_biofuel_costs.readline())

    with open(path_to_template, "r") as template_file:
        template = jinja2.Template(template_file.read())
    rendered = template.render(
        scaling_factors=scaling_factors,
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
        max_power_density=snakemake.params["max_power_density"],
        biofuel_efficiency=snakemake.params["biofuel_efficiency"],
        path_to_biofuel_costs=snakemake.input.biofuel_cost,
        path_to_result=snakemake.output[0]
    )
