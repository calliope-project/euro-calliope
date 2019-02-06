"""Creates Calliope location files."""
import jinja2
import pandas as pd


TEMPLATE = """locations:
    {{ eligibilities.index | join(', ') }}:
        techs:
            demand_elec:
            battery:
            hydrogen:
    {% for country, eligibility in eligibilities.iterrows() %}
    {{ country }}_pv_or_wind_farm:
        available_area: {{ eligibility.eligibility_onshore_wind_and_pv_km2 }}
        techs:
            open_field_pv:
                constraints:
                    resource: file=capacityfactors-open-field-pv.csv:{{ country }}
            wind_onshore:
                constraints:
                    resource: file=capacityfactors-wind-onshore.csv:{{ country }}
    {{ country }}_roof_mounted_pv:
        available_area: {{ eligibility.eligibility_rooftop_pv_km2 }}
        techs:
            roof_mounted_pv:
                constraints:
                    resource: file=capacityfactors-rooftop-pv.csv:{{ country }}
    {{ country }}_wind_offshore:
        available_area: {{ eligibility.eligibility_offshore_wind_km2 }}
        techs:
            wind_offshore:
                constraints:
                    resource: file=capacityfactors-wind-offshore.csv:{{ country }}
    {{ country }}_wind_onshore:
        available_area: {{ eligibility.eligibility_onshore_wind_km2 }}
        techs:
            wind_onshore:
                constraints:
                    resource: file=capacityfactors-wind-onshore.csv:{{ country }}
    {% endfor %}
links:
    {% for country, _ in eligibilities.iterrows() %}
    {{ country }},{{ country}}_pv_or_wind_farm:
        techs:
            free_transmission:
    {{ country }},{{ country}}_roof_mounted_pv:
        techs:
            free_transmission:
    {{ country }},{{ country}}_wind_offshore:
        techs:
            free_transmission:
    {{ country }},{{ country}}_wind_onshore:
        techs:
            free_transmission:
    {% endfor %}
"""


def construct_locations(path_to_land_eligibility_km2, path_to_result):
    """Generate a file that represents locations in Calliope."""
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        eligibilities=pd.read_csv(path_to_land_eligibility_km2, index_col=0)
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


if __name__ == "__main__":
    construct_locations(
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        path_to_result=snakemake.output[0]
    )
