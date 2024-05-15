import eurocalliopelib.utils as ec_utils
import numpy as np
import pandas as pd
import xarray as xr
from utils import jrc_idees_parser as jrc


def fill_missing_countries_years(
    eurostat_balances: pd.DataFrame,
    cat_names: pd.DataFrame,
    carrier_names: pd.DataFrame,
    jrc_subsector_demand: xr.DataArray,
) -> xr.DataArray:
    """
    For countries without relevant data in JRC_IDEES we use their
    Eurostat energy balance data to estimate future energy consumption, relative to the
    energy balance data of the 28 countries for which we do have JRC IDEES data.
    JRC data is available until 2015, so for future years data for all countries is filled
    based on Eurostat energy balances.
    Any other missing years are filled in based on average consumption of a country.
    """
    # Ensure countries follow our standard alpha3 convention.
    jrc_countries = jrc_subsector_demand["country_code"].values

    # Build eurostat annual industry balances
    eurostat_industry_balances = (
        eurostat_balances.unstack(["year", "country"])
        .groupby(
            [
                cat_names.jrc_idees.dropna().to_dict(),
                carrier_names.ind_carrier_name.dropna().to_dict(),
            ],
            level=["cat_code", "carrier_code"],
        )
        .sum(min_count=1)
        .groupby(["cat_code"])
        .sum(min_count=1)
        .stack("country")
        .rename_axis(index=["cat_name", "country_code"])
        .apply(ec_utils.tj_to_twh)
        .stack()
        .to_xarray()
    )
    # If in JRC-IDEES, then data already was adequately processed for that country
    jrc_balances = eurostat_industry_balances.sel(country_code=jrc_countries)
    # Otherwise, there will only be data for that country in annual energy balances
    nonjrc_balances = eurostat_industry_balances.drop_sel(country_code=jrc_countries)
    # Obtain share of total energy demand in missing countries per year
    nonjrc_subsector_share = nonjrc_balances / jrc_balances.sum("country_code")

    # Fill missing countries in relation to their total energy share.
    # E.g., if CHE consumes a total share of 1% -> assume it consumes 1% of total electricity
    total_annual_demand = jrc_subsector_demand.sum("country_code")
    nonjrc_subsector_demand = nonjrc_subsector_share * total_annual_demand
    all_subsector_demand = jrc_subsector_demand.combine_first(nonjrc_subsector_demand)

    # Sometimes JRC-IDEES has consumption, but Eurostat doesn't, leading to inf values
    # E.g., RO: non ferrous metals
    # Correct and fill with mean values.
    country_mean_demand = (all_subsector_demand / eurostat_industry_balances).mean(
        "year"
    )
    country_mean_demand = country_mean_demand.where(lambda x: ~np.isinf(x) & (x > 0))
    country_mean_demand = country_mean_demand.fillna(
        country_mean_demand.mean("country_code")
    )
    # Fill data where JRC says there is no consumption of any form in a country's industry subsector
    # But where the energy balances show consumption (e.g. UK, Wood and wood products)
    _to_fill = all_subsector_demand
    _filler = eurostat_industry_balances * country_mean_demand

    extra_years = [y for y in _filler.year if y > _to_fill.year.max()]
    _to_fill = xr.concat([_to_fill, _filler.sel(year=extra_years)], dim="year")
    _to_fill = _to_fill.where(_to_fill.sum("carrier_name") > 0).fillna(_filler)

    # Fill remaining missing year data
    # Interpolate middle years -> backfill older years -> forwardfill newer years
    _to_fill = _to_fill.interpolate_na(
        dim="year", use_coordinate="year", method="linear"
    )
    _to_fill = _to_fill.bfill(dim="year")
    all_filled = _to_fill.ffill(dim="year")

    # TODO: CHE has no values for "Wood and wood products" and "Transport Equipment".
    # We fill this with zeros for now (equivalent to SCEC). Should CHE data be improved/filled in another way?
    all_filled = all_filled.fillna(0)

    all_filled = jrc.standardize(all_filled, "twh")

    assert ~all_filled.isnull().any(), "Filling failed, found null values."
    assert ~np.isinf(all_filled).any(), "Filling failed, found inf values."

    return all_filled
