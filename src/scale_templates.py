"""Applies scaling factor to template files."""
import jinja2


def scale_template(path_to_template, scaling_factors, path_to_result):
    """Applies scaling factor to template files."""

    scaling_factors["specific_costs"] = scaling_factors["monetary"] / scaling_factors["power"]

    with open(path_to_template, "r") as template_file:
        template = jinja2.Template(template_file.read())
    rendered = template.render(
        scaling_factors=scaling_factors
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    scale_template(
        path_to_template=snakemake.input.template,
        scaling_factors=snakemake.params["scaling_factors"],
        path_to_result=snakemake.output[0]
    )
