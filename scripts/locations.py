"""Creates Calliope location files."""
import jinja2
import pandas as pd
import geopandas as gpd

import filters

TEMPLATE = """locations:
    {% for id, location in locations.iterrows() %}
    {{ id }}: # {{ location["name"] }}
        coordinates: {lat: {{ location.centroid.y }}, lon: {{ location.centroid.x }}}
        available_area: {{ location.eligibility_onshore_wind_and_pv_km2 * scaling_factors.area }} # {{ (1 / scaling_factors.area) | unit("km2") }} usable by onshore wind or open field pv
        techs:
            demand_elec:
            battery:
            hydrogen:
            open_field_pv:
            wind_onshore_competing:
            wind_onshore_monopoly:
                constraints:
                    energy_cap_max: {{ location.eligibility_onshore_wind_monopoly_mw * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            roof_mounted_pv:
                constraints:
                    energy_cap_max: {{ location.eligibility_rooftop_pv_mw * scaling_factors.power  }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            wind_offshore:
                constraints:
                    energy_cap_max: {{ location.eligibility_offshore_wind_mw * scaling_factors.power  }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            hydro_run_of_river:
                constraints:
                    energy_cap_max: {{ location.installed_capacity_hror_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
            hydro_reservoir:
                constraints:
                    energy_cap_max: {{ location.installed_capacity_hdam_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                    storage_cap_max: {{ location.storage_capacity_hdam_MWh * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MWh") }}
            pumped_hydro:
                constraints:
                    energy_cap_max: {{ location.installed_capacity_hphs_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                    storage_cap_max: {{ location.storage_capacity_hphs_MWh * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MWh") }}
            biofuel:
                constraints:
                    resource: {{ location.biofuel_potential_mwh_per_year / 8760 * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                    storage_cap_equals: {{ location.biofuel_potential_mwh_per_year / 2 * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MWh") }} (0.5x annual yield) # ASSUME < 1 for numerical range
    {% endfor %}
overrides:
    freeze-hydro-capacities:
        locations:
            {% for id, location in locations.iterrows() %}
            {{ id }}.techs: # {{ location["name"] }}
                hydro_run_of_river:
                    constraints:
                        energy_cap_equals: {{ location.installed_capacity_hror_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                hydro_reservoir:
                    constraints:
                        energy_cap_equals: {{ location.installed_capacity_hdam_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                        storage_cap_equals: {{ location.storage_capacity_hdam_MWh * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MWh") }}
                pumped_hydro:
                    constraints:
                        energy_cap_equals: {{ location.installed_capacity_hphs_MW * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MW") }}
                        storage_cap_equals: {{ location.storage_capacity_hphs_MWh * scaling_factors.power }} # {{ (1 / scaling_factors.power) | unit("MWh") }}
            {% endfor %}

"""


def construct_locations(path_to_shapes, path_to_land_eligibility_km2, path_to_hydro_capacities_mw,
                        path_to_biofuel_potential_mwh, flat_roof_share, maximum_installable_power_density,
                        scaling_factors, biofuel_efficiency, path_to_output_yaml, path_to_output_csv):
    """Generate a file that represents locations in Calliope."""
    locations = gpd.GeoDataFrame(
        gpd.read_file(path_to_shapes).set_index("id")
    )
    locations = (
        locations
        .assign(centroid=locations.centroid.rename("centroid"))
        .loc[:, ["name", "centroid"]]
    )
    capacities = _from_area_to_installed_capacity(
        land_eligibiligy_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        flat_roof_share=flat_roof_share,
        maximum_installable_power_density=maximum_installable_power_density
    )
    hydro_capacities = pd.read_csv(path_to_hydro_capacities_mw, index_col=0)
    biofuel = pd.read_csv(path_to_biofuel_potential_mwh, index_col=0) * biofuel_efficiency
    locations = locations.merge(
        pd.concat([capacities, hydro_capacities, biofuel], axis="columns", sort=True),
        how="left",
        left_index=True,
        right_index=True,
        validate="one_to_one"
    )
    locations = locations.assign(id=locations.index.str.replace(".", "-")).set_index("id")

    env = jinja2.Environment()
    env.filters["unit"] = filters.unit
    rendered = env.from_string(TEMPLATE).render(
        locations=locations,
        scaling_factors=scaling_factors
    )
    with open(path_to_output_yaml, "w") as result_file:
        result_file.write(rendered)
    locations.name.to_csv(path_to_output_csv, index=True, header=True)


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
        path_to_biofuel_potential_mwh=snakemake.input.biofuel,
        path_to_output_yaml=snakemake.output.yaml,
        path_to_output_csv=snakemake.output.csv,
        flat_roof_share=snakemake.params["flat_roof_share"],
        maximum_installable_power_density=snakemake.params["maximum_installable_power_density"],
        scaling_factors=snakemake.params["scaling_factors"],
        biofuel_efficiency=snakemake.params.biofuel_efficiency
    )
