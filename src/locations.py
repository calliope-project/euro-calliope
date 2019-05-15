"""Creates Calliope location files."""
import jinja2
import pandas as pd
import geopandas as gpd

FLAT_ROOF_SHARE = 0.302 # TODO add source (own publication)
MAXIMUM_INSTALLABLE_POWER_DENSITY = {
    "pv-on-tilted-roofs": 160, # [MW/km^2] from (Gagnon:2016, Klauser:2016), i.e. 16% efficiency
    "pv-on-flat-areas": 80, # [MW/km^2] from (Gagnon:2016, Klauser:2016, Wirth:2017) # FIXME also defined in renewable-techs.yaml
    "onshore-wind": 8, # [MW/km^2] from (European Environment Agency, 2009) # FIXME also defined in renewable-techs.yaml
    "offshore-wind": 15 # [MW/km^2] from (European Environment Agency, 2009)
}

TEMPLATE = """locations:
    {% for id, location in locations.iterrows() %}
    {{ id | replace(".", "-") }}:
        coordinates: {lat: {{ location.centroid.y }}, lon: {{ location.centroid.x }}}
        available_area: {{ location.eligibility_onshore_wind_and_pv_km2 * scaling_factors.area }} # [{{ 1 / scaling_factors.area }} km2] usable by onshore wind or open field pv
        techs:
            demand_elec:
            battery:
            hydrogen:
            pumped_hydro:
            hydro_run_of_river:
            open_field_pv:
            wind_onshore_competing:
            wind_onshore_monopoly:
                constraints:
                    energy_cap_max: {{ location.eligibility_onshore_wind_monopoly_mw * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MW]
            roof_mounted_pv:
                constraints:
                    energy_cap_max: {{ location.eligibility_rooftop_pv_mw * scaling_factors.power  }} # [{{ 1 / scaling_factors.power }} MW]
            wind_offshore:
                constraints:
                    energy_cap_max: {{ location.eligibility_offshore_wind_mw * scaling_factors.power  }} # [{{ 1 / scaling_factors.power }} MW]
    {% endfor %}
"""


def construct_locations(path_to_shapes, path_to_land_eligibility_km2, scaling_factors, path_to_result):
    """Generate a file that represents locations in Calliope."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id").centroid.rename("centroid")
    )
    capacities = _from_area_to_installed_capacity(
        land_eligibiligy_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        flat_roof_share=FLAT_ROOF_SHARE,
        maximum_installable_power_density=MAXIMUM_INSTALLABLE_POWER_DENSITY
    )
    locations = locations.merge(
        capacities,
        how="left",
        left_index=True,
        right_index=True,
        validate="one_to_one"
    )

    template = jinja2.Template(TEMPLATE)
    rendered = template.render(
        locations=locations,
        scaling_factors=scaling_factors
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
        path_to_shapes=snakemake.input.shapes,
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        path_to_result=snakemake.output[0],
        scaling_factors=snakemake.params["scaling_factors"]
    )
