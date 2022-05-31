import xarray as xr
import pandas as pd

from eurocalliopelib import utils

YEARS = slice(2000, 2018)


def sectoral_energy_balances(
    path_to_eurostat_energy_balances: str,
    path_to_ch_energy_balances: str,
    sector: str,
    path_to_carrier_name_mapping: str,
    category_name_mapping: dict,
    path_to_output: str
):
    """
    Blend Swiss and Eurostat data on annual energy balances for a specific subsector.
    At this stage we also merge some Eurostat subsectors together (namely "other" sectors)
    according to user-defined allocations to Euro-Calliope sectors.
    We also rename energy carriers from Eurostat codes to readable names.

    Args:
        path_to_eurostat_energy_balances (str):
            Path to xarray dataarray containing the main Eurostat energy balances.
        path_to_ch_energy_balances (str):
            Path to xarray dataarray containing Swiss-specific energy balances
            (can be sector-specific data only).
        sector (str):
            Name of Euro-Calliope "sector".
        path_to_carrier_name_mapping (str):
            Path to CSV file containing mapping from Eurostat carrier codes to readable name.
            Must contain columns mapping codes to the sector of interest
            (and, if merging with an "other" sector, the mapping for that sector).
        category_name_mapping (dict):
            If sector is not "industry":
                Keys = Euro-Calliope sector, values = list of Eurostat sector (category)
                codes over which energy balance data will be summed.
            If sector is "industry":
                Keys = Eurostat Industry subsector (category) code, values = JRC-IDEES subsector name.
        path_to_output (str):
            Path to output xarray dataarray to store the sector-specific, blended and renamed data.
    """

    eurostat_energy_balances = xr.open_dataarray(path_to_eurostat_energy_balances)
    ch_energy_balances = xr.open_dataarray(path_to_ch_energy_balances)
    all_energy_balances = utils.merge_da([eurostat_energy_balances, ch_energy_balances], "annual_energy_balances")

    carrier_name_mapping = pd.read_csv(path_to_carrier_name_mapping, index_col=0)

    if f"other-{sector}" in category_name_mapping.keys():
        # ASSUME: energy demand in "other" sectors (incl. military, fisheries, agriculture)
        # can be assigned to main Euro-Calliope "sectors" according to the type of demand
        sectoral_energy_balance = add_other_sector_to_main_sector(
            all_energy_balances, sector, f"other-{sector}", carrier_name_mapping, category_name_mapping
        )

    elif sector == "industry":
        sectoral_energy_balance = rename_carriers(
            eurostat_industry_subsectors_to_jrc_names(all_energy_balances, category_name_mapping),
            carrier_name_mapping, sector
        )

    else:
        sectoral_energy_balance = rename_carriers(
            slice_on_sector(all_energy_balances, category_name_mapping[sector]),
            carrier_name_mapping, sector
        )

    sectoral_energy_balance = sectoral_energy_balance.sel(year=YEARS)

    sectoral_energy_balance.to_netcdf(path_to_output)


def eurostat_industry_subsectors_to_jrc_names(all_energy_balances, category_name_mapping):
    return utils.rename_and_groupby(
        all_energy_balances, category_name_mapping,
        dim_name="cat_code", new_dim_name="cat_name"
    )


def slice_on_sector(all_energy_balances, category_name_mapping):
    return (
        all_energy_balances
        .sel(cat_code=category_name_mapping)
        .sum("cat_code", min_count=1, keep_attrs=True)
    )


def rename_carriers(sectoral_energy_balance, carrier_name_mapping, sector):
    sectoral_energy_balance = utils.rename_and_groupby(
        sectoral_energy_balance, carrier_name_mapping[f"{sector}_carrier_name"].dropna().to_dict(),
        dim_name="carrier_code", new_dim_name="carrier_name"
    )

    return sectoral_energy_balance


def add_other_sector_to_main_sector(
    energy_balances, main_sector, other_sector, carrier_name_mapping, category_name_mapping
):
    main_energy_balance = rename_carriers(
        slice_on_sector(energy_balances, category_name_mapping[main_sector]),
        carrier_name_mapping, main_sector
    )
    other_energy_balance = rename_carriers(
        slice_on_sector(energy_balances, category_name_mapping[other_sector]),
        carrier_name_mapping, other_sector
    ).reindex_like(main_energy_balance).fillna(0)

    summed_energy_balance = main_energy_balance + other_energy_balance

    assert (
        summed_energy_balance.groupby("carrier_name").sum(...) >=
        main_energy_balance.groupby("carrier_name").sum(...)
    ).all()

    return summed_energy_balance


if __name__ == "__main__":
    sectoral_energy_balances(
        path_to_eurostat_energy_balances=snakemake.input.eurostat_energy_balances,
        path_to_ch_energy_balances=snakemake.input.ch_energy_balances,
        sector=snakemake.wildcards.sector,
        path_to_carrier_name_mapping=snakemake.params.carrier_name_mapping,
        category_name_mapping=snakemake.params.category_name_mapping,
        path_to_output=snakemake.output[0]
    )
