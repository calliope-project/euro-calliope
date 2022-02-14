import pandas as pd

def set_capacity_limits(
    path_to_units, path_to_land_eligibility_km2, roof_shares, path_to_output, max_power_densities
):
    capacities = _from_area_to_installed_capacity(
        land_eligibility_km2=pd.read_csv(path_to_land_eligibility_km2, index_col=0),
        roof_shares=roof_shares,
        maximum_installable_power_density=max_power_densities
    )
    units = pd.read_csv(path_to_units).set_index("id")
    capacities = capacities.reindex(units.index).fillna(0)

    capacities.to_csv(path_to_output)


def _from_area_to_installed_capacity(land_eligibility_km2, roof_shares,
                                     maximum_installable_power_density):
    cap = land_eligibility_km2.copy()
    factor_rooftop = (
        maximum_installable_power_density["pv-on-flat-areas"] * roof_shares["flat"]
        + maximum_installable_power_density["pv-on-tilted-roofs"] * (1 - roof_shares["flat"])
    )
    factor_onshore = maximum_installable_power_density["onshore-wind"]
    factor_offshore = maximum_installable_power_density["offshore-wind"]
    cap["eligibility_rooftop_pv_mw"] = cap["eligibility_rooftop_pv_km2"] * factor_rooftop
    cap["eligibility_offshore_wind_mw"] = cap["eligibility_offshore_wind_km2"] * factor_offshore
    cap["eligibility_onshore_wind_monopoly_mw"] = cap["eligibility_onshore_wind_km2"] * factor_onshore
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
    set_capacity_limits(
        path_to_units=snakemake.input.units,
        path_to_land_eligibility_km2=snakemake.input.land_eligibility_km2,
        roof_shares=snakemake.params["roof_shares"],
        max_power_densities=snakemake.params["max_power_densities"],
        path_to_output=snakemake.output[0]
    )
