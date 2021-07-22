"""Create override allowing load shedding."""
import jinja2
import pandas as pd

TEMPLATE = """overrides:
    load-shedding:
        locations:
            {% for id, location in locations.iterrows() %}
            {{ id | replace(".", "-") }}.techs.load_shedding:
            {% endfor %}
"""


def load_shedding(path_to_units, path_to_result):
    """Generate a file that allows load shedding."""
    locations = pd.read_csv(path_to_units, index_col=0)
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        locations=locations
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    load_shedding(
        path_to_units=snakemake.input.units,
        path_to_result=snakemake.output[0]
    )
