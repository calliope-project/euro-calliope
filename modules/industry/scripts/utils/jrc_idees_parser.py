from typing import Optional, Union

import numpy as np
import xarray as xr

STANDARD_COORDS = ["cat_name", "year", "country_code", "carrier_name"]


def check_units(jrc_energy: xr.Dataset, jrc_prod: xr.DataArray) -> None:
    """Check that the JRC data is in the right units."""
    for var in jrc_energy:
        assert jrc_energy[var].attrs["units"].lower() == "twh"
    assert jrc_prod.attrs["units"].lower() == "kt"


def standardize(
    da: xr.DataArray, units: str, name: Optional[str] = None
) -> xr.DataArray:
    """Ensure JRC processing standard is met.

    Three requirements:
    1. coordinates must follow a given naming convention.
    2. attribute for units must be set.
    3. name must be set.

    Args:
        da (xr.DataArray): target xarray object to standardize.
        units (str): units to set as attribute.
        name (Optional[str], optional): name to set. Defaults to None.

    Raises:
        KeyError: If standard coordinates cannot be set.
        ValueError: If name cannot be set.

    Returns:
        xr.DataArray: standardized xarray object.
    """
    # 1. coordinate standard
    removed_coords = set(da.coords).difference(STANDARD_COORDS)
    removed_dims = removed_coords.intersection(da.dims)
    if removed_dims:
        raise KeyError(f"Cannot ensure standard coordinates for {removed_dims}.")
    da = da.drop(removed_coords)
    # 2. units standard
    da = da.assign_attrs(units=units)
    # 3. name standard
    if da.name is None and name is None:
        raise ValueError("Name must be set for DataArrays without it.")
    if name:
        da.name = name

    return da


def get_section_final_intensity_auxiliary_electric(
    section: str,
    material: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.DataArray,
    fill_empty: bool = False,
) -> xr.DataArray:
    """Wrapper for auxiliary electrical processes."""
    auxiliaries = ["Lighting", "Air compressors", "Motor drives", "Fans and pumps"]
    auxiliary_intensity = sum(
        get_section_subsection_final_intensity(
            section, aux, material, "Electricity", jrc_energy, jrc_prod, fill_empty
        )
        for aux in auxiliaries
    )

    return auxiliary_intensity


def get_section_subsection_final_intensity(
    section: str,
    subsection: str,
    material: str,
    carrier_name: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.DataArray,
    fill_empty: bool = False,
) -> xr.DataArray:
    """Get final energy intensity of a given JRC section/subsection/material."""
    # Extract relevant section and subsection data.
    final_demand = jrc_energy["final"].sel(section=section, subsection=subsection)
    useful_demand = jrc_energy["useful"].sel(section=section, subsection=subsection)
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
    final_intensity = standardize(final_intensity, "twh/kt", name="final_intensity")
    final_intensity.name = "final_intensity"

    assert final_intensity.sum() >= useful_intensity.sum(), "Creating energy!"

    return final_intensity.fillna(0)


def get_section_subsection_useful_intensity(
    section: str,
    subsection: str,
    material: str,
    jrc_energy: xr.Dataset,
    jrc_prod: xr.DataArray,
) -> xr.DataArray:
    """Get useful energy intensity of a given section/subsection/material."""
    useful_demand = jrc_energy["useful"].sel(section=section, subsection=subsection)
    production = jrc_prod.sel(produced_material=material)
    useful_intensity = (useful_demand / production).sum("carrier_name")

    # Prettify
    useful_intensity = standardize(useful_intensity, "twh/kt", name="useful_intensity")

    assert ~np.isinf(
        useful_intensity
    ).any(), f"Zero division ocurred for {section}.{subsection} and {material}!"

    return useful_intensity.fillna(0)


def replace_carrier_final_demand(
    carrier: str,
    jrc_energy: xr.Dataset,
) -> xr.DataArray:
    """Assume that all demand that could consume a carrier is completely met by that carrier.

    Uses the end-use demand and conversion efficiencies as proxy for filling.

    Args:
        carrier (str): carrier to use as replacement.
        jrc_energy (xr.Dataset): JRC energy dataset.

    Returns:
        xr.DataArray: final demand data met using the specified carrier.
    """
    useful_dem_tot = jrc_energy["useful"].sum(
        "carrier_name"
    )  # MC note: including all carriers

    carrier_tot = jrc_energy.sel(
        carrier_name=carrier
    )  # MC note: only carrier_ralted things here
    carrier_eff = carrier_tot["useful"] / carrier_tot["final"]

    # ASSUME: fill NaNs (where there is demand, but no consumption in that country)
    # First by country avg. (all years), then by year avg. (all countries).
    carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="year"))
    carrier_eff = carrier_eff.fillna(carrier_eff.mean(dim="country_code"))

    carrier_final_demand = (
        useful_dem_tot / carrier_eff
    )  # MC note: this basically means all carriers demand added together divided by specific carrier efficiency? But this formula only works for the entries that relate to the defined carrier
    carrier_final_demand = carrier_final_demand.dropna(
        dim="section", how="all"
    )  # MC note: added to exclude nans effectively
    carrier_final_demand = carrier_final_demand.dropna(
        dim="subsection", how="all"
    )  # MC note: this means that the subsectors demand that cannot be met by the defined carrier are removed. For chemicals industry, subsection shrinks from 78 to 9

    # Prettify
    carrier_final_demand = carrier_final_demand.drop(["carrier_name"])
    carrier_final_demand = carrier_final_demand.assign_attrs(units="twh")
    carrier_final_demand.name = "final"
    # assert useful_dem_tot.sum() < carrier_final_demand.sum(), "Creating energy!" # MC note: but useful_dem_tot includes all carrier demands while carrier_final_demand only contains the demand related to one specific carrier. Surely this is fair comparison?
    # MC note: commented out because this is not correct for chemicals industry. but need to coordinate with Ivan regarding how this function should be written.
    return carrier_final_demand


def convert_subsection_demand_to_carrier(
    jrc_energy: xr.Dataset,
    subsection: Union[str, list[str]],
    demand_type: str = "useful",
) -> xr.DataArray:
    """Converts a subsection into a carrier by aggregating all demand for it.

    Args:
        jrc_energy (xr.Dataset): JRC energy dataset.
        subsection (Union[str, list[str]]): subsections to transform into carriers.
        demand_type (str, optional): demand type of the carrier. Defaults to "useful".

    Returns:
        xr.DataArray: new carrier data in standard coordinate format.
    """
    subsec_useful_dem = jrc_energy.sel(subsection=subsection)[demand_type]

    new_carrier_useful_dem = subsec_useful_dem.sum(["carrier_name", "section"])
    new_carrier_useful_dem = new_carrier_useful_dem.rename({
        "subsection": "carrier_name"
    })

    return standardize(new_carrier_useful_dem, "twh")
