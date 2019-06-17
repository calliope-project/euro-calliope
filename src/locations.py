"""Creates Calliope location files."""
import jinja2
import pandas as pd
import geopandas as gpd


TEMPLATE = """locations:
    {% for id, location in locations.iterrows() %}
    {{ id | replace(".", "-") }}:
        coordinates: {lat: {{ location.centroid.y }}, lon: {{ location.centroid.x }}}
        available_area: {{ location.eligibility_onshore_wind_and_pv_km2 * scaling_factors.area }} # [{{ 1 / scaling_factors.area }} km2] usable by onshore wind or open field pv
        techs:
            demand_elec:
            battery:
            hydrogen:
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
            hydro_run_of_river:
                constraints:
                    energy_cap_equals: {{ location.installed_capacity_hror_MW * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MW]
            hydro_reservoir:
                constraints:
                    energy_cap_equals: {{ location.installed_capacity_hdam_MW * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MW]
                    storage_cap_equals: {{ location.storage_capacity_hdam_MWh * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MWh]
            pumped_hydro:
                constraints:
                    energy_cap_equals: {{ location.installed_capacity_hphs_MW * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MW]
                    storage_cap_equals: {{ location.storage_capacity_hphs_MWh * scaling_factors.power }} # [{{ 1 / scaling_factors.power }} MWh]
    {% endfor %}
"""


def construct_locations(path_to_shapes, path_to_land_eligibility_km2, path_to_hydro_capacities_mw,
                        flat_roof_share, maximum_installable_power_density, scaling_factors, path_to_result):
    """Generate a file that represents locations in Calliope."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id").centroid.rename("centroid")
    )
    capacities = _from_area_to_installed_capacity(
        land_eligibiligy_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        flat_roof_share=flat_roof_share,
        maximum_installable_power_density=maximum_installable_power_density
    )
    hydro_capacities = pd.read_csv(path_to_hydro_capacities_mw, index_col=0)
    locations = locations.merge(
        pd.concat([capacities, hydro_capacities], axis="columns"),
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
        path_to_hydro_capacities_mw=snakemake.input.hydro_capacities,
        path_to_result=snakemake.output[0],
        flat_roof_share=snakemake.params["flat_roof_share"],
        maximum_installable_power_density=snakemake.params["maximum_installable_power_density"],
        scaling_factors=snakemake.params["scaling_factors"]
    )
