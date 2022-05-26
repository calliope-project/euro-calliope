import pandas as pd
import xarray as xr
import numpy as np

xr.set_options(keep_attrs=True)

def annual_heat_demand(
    path_to_end_use_energy_supply, path_to_annual_energy_balances,
    heat_tech_efficiency_overrides, path_to_output
):
    end_use_energy_supply = xr.open_dataarray(path_to_end_use_energy_supply)
    annual_energy_balances = xr.open_dataarray(path_to_annual_energy_balances)

    end_use_energy_supply = select_heat_end_uses(end_use_energy_supply)
    end_use_energy_supply = (
        end_use_energy_supply
        .reindex(end_use=["cooking", "water_heat", "space_heat"])
        .dropna("end_use", how="all")
    )

    efficiencies = pd.DataFrame(
        data=1,
        index=end_use_energy_supply.carrier_name.to_index(),
        columns=end_use_energy_supply.end_use.to_index()

    )
    efficiency_overrides = pd.DataFrame(heat_tech_efficiency_overrides)
    efficiency_overrides.columns = efficiency_overrides.columns.str.replace("-", "_")

    efficiencies.update(efficiency_overrides)
    efficiencies_da = efficiencies.stack().to_xarray()

    end_use_demand = (end_use_energy_supply * efficiencies_da).sum("carrier_name", min_count=1)

    end_use_demand = fill_missing_data_in_years(
        end_use_demand, annual_energy_balances, fill_method="mean"
    )

    # test: demand is at least the total energy consumption multiplied by the worst efficiency

    min_demand = (end_use_energy_supply * efficiencies_da.min("carrier_name")).sum("carrier_name")
    assert (end_use_demand >= min_demand).all()

    end_use_demand.to_netcdf(path_to_output)


def current_annual_electricity_consumption(
    path_to_end_use_energy_supply, path_to_annual_energy_balances, path_to_output
):
    """
    Get annual energy consumption coming from electricity, to remove later from the
    electricity demand profile.
    Gaps in electricity consumption are filled before saving, based on annual energy
    consumption
    """
    end_use_energy_supply = xr.open_dataarray(path_to_end_use_energy_supply)
    annual_energy_balances = xr.open_dataarray(path_to_annual_energy_balances)

    end_use_energy_supply = select_heat_end_uses(end_use_energy_supply)
    electricity_consumption = (
        end_use_energy_supply
        .drop_sel(end_use="end_use_electricity", errors="ignore")
        .reindex(carrier_name=["electricity", "direct_electric", "heat_pump"])
        .sum("carrier_name", min_count=1)
    )

    electricity_consumption = fill_missing_data_in_years(
        electricity_consumption, annual_energy_balances, fill_method="first"
    )

    electricity_consumption.to_netcdf(path_to_output)


def select_heat_end_uses(da):
    return (
        da
        .reindex(end_use=["cooking", "water_heat", "space_heat"])
        .dropna("end_use", how="all")
    )


def fill_missing_data_in_years(end_use_demand, annual_energy_balance, fill_method):
    # ASSUME: A sector's total annual energy consumption per country and year
    # scales with that sector's building demand.
    # We use this assumption to fill years with missing demand data.

    # get total energy supply from Eurostat annual energy balances
    country_energy_balance = annual_energy_balance.sum("carrier_name").where(lambda x: x > 0)

    # get end use demand per unit total energy supply, averaged over all years
    average_annual_country_supply_fraction = (end_use_demand / country_energy_balance).mean("year")

    # scale total energy supply by average end use demand per unit to get a proxy for
    # end use demand in any years without it.
    end_use_demand = end_use_demand.fillna(
        average_annual_country_supply_fraction * country_energy_balance
    )

    # ASSUME: as a last resort, filling demand data directly from other years is
    # feasible, according to the given fill method, e.g. mean/first/last/max/min

    end_use_demand = end_use_demand.fillna(
        getattr(end_use_demand.groupby("country_code"), fill_method)(...)
    )

    return end_use_demand


if __name__ == "__main__":
    if snakemake.wildcards.demand_type == "heat-demand":
        annual_heat_demand(
            path_to_end_use_energy_supply=snakemake.input.end_use_energy_supply,
            path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
            heat_tech_efficiency_overrides=snakemake.params.heat_tech_efficiency_overrides,
            path_to_output=snakemake.output[0]
        )
    elif snakemake.wildcards.demand_type == "current-electricity-consumption":
        current_annual_electricity_consumption(
            path_to_end_use_energy_supply=snakemake.input.end_use_energy_supply,
            path_to_annual_energy_balances=snakemake.input.annual_energy_balances,
            path_to_output=snakemake.output[0]
        )
