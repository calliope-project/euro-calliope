"""Creates Calliope location files."""
import jinja2
import pandas as pd

FLAT_ROOF_SHARE = 0.302 # TODO add source (own publication)
MAXIMUM_INSTALLABLE_POWER_DENSITY = {
    "pv-on-tilted-roofs": 160, # [MW/km^2] from (Gagnon:2016, Klauser:2016), i.e. 16% efficiency
    "pv-on-flat-areas": 80, # [MW/km^2] from (Gagnon:2016, Klauser:2016, Wirth:2017) # FIXME also defined in renewable-techs.yaml
    "onshore-wind": 8, # [MW/km^2] from (European Environment Agency, 2009) # FIXME also defined in renewable-techs.yaml
    "offshore-wind": 15 # [MW/km^2] from (European Environment Agency, 2009)
}

TEMPLATE = """locations:
    {{ eligibilities.index | join(', ') }}:
        techs:
            demand_elec:
            battery:
            hydrogen:
            open_field_pv:
            wind_onshore_competing:
    {% for country, eligibility in eligibilities.iterrows() %}
    {{ country }}:
        available_area: {{ eligibility.eligibility_onshore_wind_and_pv_km2 }} # [km2] usable by onshore wind or open field pv
        techs:
            wind_onshore_monopoly:
                constraints:
                    energy_cap_max: {{ eligibility.eligibility_onshore_wind_monopoly_mw }} # [MW]
            roof_mounted_pv:
                constraints:
                    energy_cap_max: {{ eligibility.eligibility_rooftop_pv_mw }} # [MW]
            wind_offshore:
                constraints:
                    energy_cap_max: {{ eligibility.eligibility_offshore_wind_mw }} # [MW]
    {% endfor %}
"""


def construct_locations(path_to_land_eligibility_km2, path_to_result):
    """Generate a file that represents locations in Calliope."""
    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        eligibilities=_from_area_to_installed_capacity(
            land_eligibiligy_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
            flat_roof_share=FLAT_ROOF_SHARE,
            maximum_installable_power_density=MAXIMUM_INSTALLABLE_POWER_DENSITY
        )
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


def _from_area_to_installed_capacity(land_eligibiligy_km2, flat_roof_share,
                                     maximum_installable_power_density):
    cap = land_eligibiligy_km2.copy()
    factor_rooftop = (
        maximum_installable_power_density["pv-on-flat-areas"] * flat_roof_share
        + maximum_installable_power_density["pv-on-tilted-roofs"] * (1 - flat_roof_share)
    )
    factor_onshore = maximum_installable_power_density["onshore-wind"]
    factor_offshore = maximum_installable_power_density["offshore-wind"]
    cap["eligibility_rooftop_pv_mw"] = cap["eligibility_rooftop_pv_km2"] * factor_rooftop
    cap["eligibility_offshore_wind_mw"] = cap["eligibility_offshore_wind_km2"] * factor_offshore
    cap["eligibility_onshore_wind_monopoly_mw"] = cap["eligibility_onshore_wind_km2"] * factor_onshore
    return cap


if __name__ == "__main__":
    construct_locations(
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        path_to_result=snakemake.output[0]
    )
