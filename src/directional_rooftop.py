"""Creates Calliope location files."""
import jinja2
import pandas as pd
import geopandas as gpd

import filters

TEMPLATE = """overrides:
    directional-rooftop-pv:
        techs:
            roof_mounted_pv_s_flat:
                essentials:
                    name: Roof mounted PV south-facing and flat
                    parent: pv_on_roof
                constraints:
                    resource: file=capacityfactors-rooftop-pv-s-flat.csv
                    resource_unit: energy_per_cap
            roof_mounted_pv_n:
                essentials:
                    name: Roof mounted PV north-facing
                    parent: pv_on_roof
                constraints:
                    resource: file=capacityfactors-rooftop-pv-n.csv
                    resource_unit: energy_per_cap
            roof_mounted_pv_e_w:
                essentials:
                    name: Roof mounted PV east-facing and west-facing
                    parent: pv_on_roof
                constraints:
                    resource: file=capacityfactors-rooftop-pv-e-w.csv
                    resource_unit: energy_per_cap
        locations:
            {% for id, location in locations.iterrows() %}
            {{ id | replace(".", "-") }}:
                techs:
                    roof_mounted_pv:
                        exists: False
                    roof_mounted_pv_s_flat:
                        constraints:
                            energy_cap_max: {{ location.eligibility_rooftop_pv_s_flat_mw * scaling_factors.power  }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                    roof_mounted_pv_n:
                        constraints:
                            energy_cap_max: {{ location.eligibility_rooftop_pv_n_mw * scaling_factors.power  }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                    roof_mounted_pv_e_w:
                        constraints:
                            energy_cap_max: {{ location.eligibility_rooftop_pv_e_w_mw * scaling_factors.power  }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            {% endfor %}
"""


def directional_rooftop(path_to_shapes, path_to_land_eligibility_km2,
                        roof_shares, maximum_installable_power_density,
                        scaling_factors, path_to_result):
    """Generate a file that represents locations in Calliope."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id").centroid.rename("centroid")
    )
    capacities = _from_area_to_installed_capacity(
        land_eligibiligy_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        roof_shares=roof_shares,
        maximum_installable_power_density=maximum_installable_power_density
    )
    locations = locations.merge(
        pd.concat([capacities], axis="columns"),
        how="left",
        left_index=True,
        right_index=True,
        validate="one_to_one"
    )
    env = jinja2.Environment()
    env.filters["unit"] = filters.unit
    rendered = env.from_string(TEMPLATE).render(
        locations=locations,
        scaling_factors=scaling_factors
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(rendered)


def _from_area_to_installed_capacity(land_eligibiligy_km2, roof_shares,
                                     maximum_installable_power_density):
    cap = land_eligibiligy_km2.copy()
    cap["eligibility_rooftop_pv_n_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["N"]
        * maximum_installable_power_density["pv-on-tilted-roofs"]
    )
    cap["eligibility_rooftop_pv_e_w_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * (roof_shares["E"] + roof_shares["W"])
        * maximum_installable_power_density["pv-on-tilted-roofs"]
    )
    cap["eligibility_rooftop_pv_s_flat_mw"] = (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["S"]
        * maximum_installable_power_density["pv-on-tilted-roofs"]
    ) + (
        cap["eligibility_rooftop_pv_km2"]
        * roof_shares["flat"]
        * maximum_installable_power_density["pv-on-flat-areas"]
    )
    return cap


if __name__ == "__main__":
    directional_rooftop(
        path_to_shapes=snakemake.input.shapes,
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        path_to_result=snakemake.output[0],
        roof_shares=snakemake.params["roof_shares"],
        maximum_installable_power_density=snakemake.params["maximum_installable_power_density"],
        scaling_factors=snakemake.params["scaling_factors"],
    )
