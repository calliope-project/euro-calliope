import xarray as xr
import numpy as np

from eurocalliopelib import utils

GAP_FILLING_COUNTRIES = {
    "balkans": {
        "countries_with_data": ["BGR", "HRV", "HUN", "ROU", "GRC"],
        "countries_without_data": ["ALB", "BIH", "MNE", "MKD", "SRB"]
    },
    "nordics": {
        "countries_with_data": ["SWE", "FIN", "DNK"],
        "countries_without_data": ["NOR", "ISL"]
    },
    "alpine": {
        "countries_with_data": ["AUT", "ITA", "DEU", "FRA"],
        "countries_without_data": ["CHE"]
    }
}

JRC_INDUSTRY_CARRIER_NAMES = {
    "Biomass": "biofuel",
    "Coke": "solid_fossil",
    "Derived gases": "gas",
    "Diesel oil": "oil",
    "Diesel oil (incl. biofuels)": "oil",
    "Electricity": "electricity",
    "LPG": "oil",
    "Naphtha": "oil",
    "Natural gas": "gas",
    "Natural gas (incl. biogas)": "gas",
    "Other liquids": "oil",
    "Refinery gas": "gas",
    "Residual fuel oil": "oil",
    "Solar and geothermal": "renewable_heat",
    "Solids": "solid_fossil",
    "Steam distributed": "heat"
}


def commercial_energy_supply(
    path_to_ch_building_heat_energy_supply: str,
    path_to_annual_energy_balances: str,
    path_to_jrc_idees_sector_energy_supply: str,
    countries: str,
    path_to_output: str
):
    """
    Combine allocation of supply to industry building consumption end uses
    from JRC-IDEES with total supply of energy to commercial buildings from Eurostat.
    Also add Swiss allocations (already processed in another script).

    Args:
        path_to_ch_building_heat_energy_supply (str):
            Path to xarray Dataarray containing pre-processed Swiss allocation of fuels
            to commercial building end uses.
        path_to_annual_energy_balances (str):
            Path to xarray Dataarray containing Eurostat annual energy balances.
        path_to_jrc_idees_sector_energy_supply (str):
            Path to xarray Dataarray containing JRC-IDEES commercial (a.k.a. tertiary)
            sector end-use consumption and demand.
        path_to_output (str):
            Path to output xarray Dataarray.
    """

    # Get annual energy balance data for the commercial sector
    energy_balance = xr.open_dataarray(path_to_annual_energy_balances)

    jrc_end_use_energy_supply = (
        xr.open_dataarray(path_to_jrc_idees_sector_energy_supply)
        .sel(energy="consumption")
    )
    country_codes = [utils.convert_country_code(country) for country in countries]
    # Map JRC end uses to annual commercial demand
    mapped_end_uses = map_jrc_to_eurostat(energy_balance, jrc_end_use_energy_supply)
    # We run it through the wash again to fill in any gaps from freshly added countries
    mapped_end_uses = map_jrc_to_eurostat(energy_balance, mapped_end_uses)

    # ASSUME: all ambient heat in Eurostat is consumed by heat pumps in space heating.
    # JRC data only refers to heat pumps for heating in space heating, ommitting the
    # consumption of ambient heat. Since we use ambient heat at a later point, we
    # add it directly from Eurostat here.
    mapped_end_uses = utils.merge_da([
        mapped_end_uses,
        energy_balance.sel(carrier_name=["ambient_heat"]).expand_dims(end_use=["space_heat"])
    ])

    check_no_data_loss(mapped_end_uses, energy_balance.sel(country_code=country_codes))

    ch_household_end_use_energy_balance = xr.open_dataarray(path_to_ch_building_heat_energy_supply)

    # Override gap-filled swiss data with data directly from swiss statistics
    commercial_end_use_energy_balance_incl_ch = utils.merge_da([
        mapped_end_uses.drop_sel(country_code="CHE").expand_dims(cat_name=["commercial"]),
        ch_household_end_use_energy_balance.sel(cat_name=["commercial"])
    ])

    commercial_end_use_energy_balance_incl_ch.sel(country_code=country_codes).to_netcdf(path_to_output)


def industry_energy_supply(
    path_to_annual_energy_balances: str,
    path_to_jrc_idees_sector_energy_supply: str,
    countries: list,
    path_to_output: str
):
    """
    Combine allocation of supply to industry building consumption end uses
    from JRC-IDEES with total supply of energy to industrial buildings from Eurostat.

    Args:
        path_to_annual_energy_balances (str):
            Path to xarray Dataarray containing industry subsector components of
            the Eurostat annual energy balances.
        path_to_jrc_idees_sector_energy_supply (str):
            Path to xarray Dataarray containing JRC-IDEES industry subsector
            end-use consumption and demand.
        path_to_output (str):
            Path to output xarray Dataarray.
    """

    # Get annual energy balance data for industry subsectors
    energy_balance = xr.open_dataarray(path_to_annual_energy_balances)

    jrc_end_use_energy_supply = (
        xr.open_dataarray(path_to_jrc_idees_sector_energy_supply)
        .sel(energy="consumption")
    )

    country_codes = [utils.convert_country_code(country) for country in countries]

    # We don"t care about subsubsectors , so we sum over them
    jrc_end_use_energy_supply = jrc_end_use_energy_supply.sum("section")
    jrc_end_use_energy_supply = utils.rename_and_groupby(
        jrc_end_use_energy_supply, JRC_INDUSTRY_CARRIER_NAMES, dim_name="carrier_name"
    )
    # Separate energy into space heat and *not* space heat.
    # ASSUME: all industry low enthalpy heat demand is for space heat.
    jrc_end_use_energy_supply = utils.merge_da([
        jrc_end_use_energy_supply.sel(subsection="Low enthalpy heat").drop("subsection").expand_dims(end_use=["space_heat"]),
        jrc_end_use_energy_supply.drop_sel(subsection="Low enthalpy heat").sum("subsection").expand_dims(end_use=["not_space_heat"])
    ])

    # Map JRC end uses to annual industry demand.
    mapped_end_uses = map_jrc_to_eurostat(energy_balance, jrc_end_use_energy_supply)
    # We run it through the wash again to fill in any gaps from freshly added countries
    mapped_end_uses = map_jrc_to_eurostat(energy_balance, mapped_end_uses)

    # ASSUME: all ambient heat in Eurostat is consumed by heat pumps in space heating.
    # JRC data only refers to heat pumps for heating in space heating, ommitting the
    # consumption of ambient heat. Since we use ambient heat at a later point, we
    # add it directly from Eurostat here.
    industry_end_use_energy_balance = utils.merge_da([
        mapped_end_uses,
        energy_balance.sel(carrier_name=["ambient_heat"]).expand_dims(end_use=["space_heat"])
    ])
    check_no_data_loss(industry_end_use_energy_balance, energy_balance.sel(country_code=country_codes))

    industry_end_use_energy_balance.sel(country_code=country_codes).to_netcdf(path_to_output)


def map_jrc_to_eurostat(energy_balance, jrc_end_use):
    """
    JRC-IDEES provides us with the distribution of supplied fuels to end uses in different sectors.
    These end uses include space heating, water heating, and cooking.
    We take the fraction of supplied fuel to each end use and apply it to the total
    fuel supplied to the a sector according to Eurostat.

    ASSUME: Total energy balances from Eurostat are more accurate than those from JRC-IDEES.

    Args:
        energy_balance (xr.DataArray): pre-processed sectoral Eurostat energy balances
        jrc_end_use (xr.DataArray): pre-processed sectoral JRC-IDEES end use energy consumption

    Returns:
        xr.DataArray: Sector fuel consumption by end use.
    """

    # ASSUME: missing end use allocation of an energy carrier in JRC-IDEES is assigned
    # to space heating
    missed_energy = (
        jrc_end_use.sum("end_use", min_count=1).where(lambda x: x > 0).isnull() &
        energy_balance.where(lambda x: x > 0).notnull()
    )
    jrc_end_use.loc[{"end_use": "space_heat"}] = (
        jrc_end_use.sel(end_use="space_heat")
        .where(lambda x: x > 0)
        .fillna(energy_balance.where(missed_energy))
    )

    # Get percentage contribution of each end use (cooking, water heating, space heating, etc.)
    # to total final energy consumption in the given sector
    jrc_end_use_fraction = (
        jrc_end_use / jrc_end_use.sum("end_use", min_count=1)
    )
    # ASSUME: Missing years (namely 2016-2018) have the same percentage contribution
    # as the average across all other years
    missing_years = energy_balance.year.to_index().difference(jrc_end_use.year.to_index())
    for yr in missing_years:
        missing_year_filler = jrc_end_use_fraction.mean("year").expand_dims(year=[yr])
        missing_year_filler = missing_year_filler / missing_year_filler.sum("end_use")
        jrc_end_use_fraction = utils.merge_da([jrc_end_use_fraction, missing_year_filler])

    # ASSUME: missing countries have the same percentage contribution as their neighbours
    for country_group_mappings in GAP_FILLING_COUNTRIES.values():
        country_averages = (
            jrc_end_use_fraction
            .sel(country_code=country_group_mappings["countries_with_data"])
            .mean("country_code")
        )
        country_averages = country_averages / country_averages.sum("end_use")
        for country in country_group_mappings["countries_without_data"]:
            if country not in jrc_end_use_fraction.country_code:
                jrc_end_use_fraction = utils.merge_da(
                    [jrc_end_use_fraction, country_averages.expand_dims(country_code=[country])]
                )
    # multiply eurostat energy balances by the JRC-IDEES percentage contributions
    mapped_end_uses = jrc_end_use_fraction * energy_balance

    return mapped_end_uses


def check_no_data_loss(mapped_end_uses, original_energy_balance):
    # check we haven"t lost any data by mapping end uses

    assert np.allclose(
        mapped_end_uses
        .sum("end_use", min_count=1)
        .where(lambda x: x > 0)
        .transpose(*original_energy_balance.dims)
        .reindex_like(original_energy_balance),
        original_energy_balance.where(lambda x: x > 0),
        equal_nan=True
    )


if __name__ == "__main__":
    if snakemake.wildcards.building_sector == "commercial":
        commercial_energy_supply(
            path_to_ch_building_heat_energy_supply=snakemake.input.ch_building_heat_energy_supply,
            path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
            path_to_jrc_idees_sector_energy_supply=snakemake.input.jrc_idees_sector_energy_supply,
            countries=snakemake.params.countries,
            path_to_output=snakemake.output[0],
        )
    if snakemake.wildcards.building_sector == "industry":
        industry_energy_supply(
            path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
            path_to_jrc_idees_sector_energy_supply=snakemake.input.jrc_idees_sector_energy_supply,
            countries=snakemake.params.countries,
            path_to_output=snakemake.output[0],
        )
