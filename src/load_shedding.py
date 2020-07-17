"""Create override allowing load shedding."""
import jinja2
import geopandas as gpd

TEMPLATE = """overrides:
    load-shedding:
        locations:
            {% for id, location in locations.iterrows() %}
            {{ id | replace(".", "-") }}.techs.load_shedding:
            {% endfor %}
"""


def load_shedding(path_to_shapes, path_to_result):
    """Generate a file that allows load shedding."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id")
    )
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        locations=locations
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    load_shedding(
        path_to_shapes=snakemake.input.shapes,
        path_to_result=snakemake.output[0]
    )
