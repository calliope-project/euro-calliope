import xarray as xr

STANDARD_COORDS = ["cat_name", "year", "country_code", "carrier_name"]


def ensure_standard_coordinates(ds: xr.Dataset):
    """Remove all coordinates that do not match a predefined standard."""
    removed_coords = [i for i in ds.coords if i not in STANDARD_COORDS]
    removed_dims = [i for i in ds.dims if i in removed_coords]
    if any(removed_dims):
        raise ValueError(f"Cannot ensure standard coordinates for {removed_dims}.")

    ds = ds.drop(removed_coords)
    return ds


def get_auxiliary_electric_final_intensity(
    section: str,
    material: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.Dataset,
    fill_empty: bool = False,
) -> xr.Dataset:
    """Wrapper for auxiliary electrical processes."""
    auxiliaries = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    auxiliary_intensity = sum(
        get_subsection_final_intensity(
            section, aux, material, "Electricity", jrc_energy, jrc_prod, fill_empty
        )
        for aux in auxiliaries
    )

    return auxiliary_intensity


def get_subsection_final_intensity(
    section: str,
    subsection: str,
    material: str,
    carrier_name: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.Dataset,
    fill_empty: bool = False,
) -> xr.Dataset:
    """Get final energy intensity of a given JRC section/subsection/material."""
    # Extract relevant section and subsection data.
    final_demand = jrc_energy.sel(
        energy="consumption", section=section, subsection=subsection
    )
    useful_demand = jrc_energy.sel(
        energy="demand", section=section, subsection=subsection
    )
    production = jrc_prod.sel(produced_material=material)

    total_eff = useful_demand / final_demand
    carrier_eff = total_eff.where(total_eff > 0).sel(carrier_name=carrier_name)
    if fill_empty:
        # First by country avg. (all years), then by year avg. (all countries).
        carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="year"))
        carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="country_code"))

    # Get the useful energy intensity of all production (e.g., twh/kt_steel)
    useful_intensity = useful_demand.sum(dim="carrier_name") / production
    # Then reconstruct final intensity.
    final_intensity = useful_intensity / carrier_eff

    # Prettify
    final_intensity = ensure_standard_coordinates(final_intensity)
    final_intensity["value"].attrs["units"] = "twh/kt"

    assert final_intensity >= useful_intensity, "Creating energy!"

    return final_intensity.fillna(0)


def get_subsection_useful_intensity(
    section: str,
    subsection: str,
    material: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.Dataset,
):
    """Get useful energy intensity of a given section/subsection/material."""
    useful_demand = jrc_energy.sel(
        energy="demand", section=section, subsection=subsection
    )
    production = jrc_prod.sel(produced_material=material)
    useful_intensity = (useful_demand / production).sum("carrier_name")

    # Prettify
    useful_intensity = ensure_standard_coordinates(useful_intensity)
    useful_intensity["value"].attrs["units"] = "twh/kt"

    return useful_intensity.fillna(0)


def get_carrier_final_demand(
    carrier: str,
    all_useful_demand: xr.Dataset,
    jrc_energy: xr.Dataset,
    fill_empty: bool = True,
) -> xr.Dataset:
    """
    Get demand for a specific carrier, assuming all end use demand that could consume
    that carrier are completely met by that carrier.
    """
    carrier_total = jrc_energy.sel(carrier_name=carrier)
    carrier_eff = carrier_total.sel(energy="demand") / carrier_total.sel(
        energy="consumption"
    )
    # Fill NaNs (where there is demand, but no consumption in that country)
    if fill_empty:
        # First by country avg. (all years), then by year avg. (all countries).
        carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="year"))
        carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="country_code"))

    carrier_final_demand = all_useful_demand / carrier_eff
    carrier_final_demand = carrier_final_demand.dropna(dim="subsection", how="all")
    carrier_final_demand = carrier_final_demand.drop(["energy", "carrier_name"])

    return carrier_final_demand
