"""Applies config parameters to template files."""
import jinja2


def parameterise_template(path_to_template, scaling_factors, max_power_density, path_to_result):
    """Applies config parameters to template files."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    with open(path_to_template, "r") as template_file:
        template = jinja2.Template(template_file.read())
    rendered = template.render(
        scaling_factors=scaling_factors,
        max_power_density=max_power_density
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    parameterise_template(
        path_to_template=snakemake.input.template,
        scaling_factors=snakemake.params["scaling_factors"],
        max_power_density=snakemake.params["max_power_density"],
        path_to_result=snakemake.output[0]
    )
