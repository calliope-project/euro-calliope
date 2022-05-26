import xarray as xr
import pandas as pd

from eurocalliopelib import utils

EUROSTAT_SECTOR_MAPPING = {  # industry subsector name mappings come from `path_to_category_names`
    "household": ["FC_OTH_HH_E"],
    "commercial": ["FC_OTH_CP_E"],
    "other": ["FC_OTH_AF_E", "FC_OTH_FISH_E", "FC_OTH_NSP_E"]
}
YEARS = slice(2000, 2018)


def sectoral_energy_balances(
    path_to_eurostat_energy_balances, path_to_ch_energy_balances, sector,
    path_to_carrier_names, path_to_category_names, path_to_output
):
    eurostat_energy_balances = xr.open_dataarray(path_to_eurostat_energy_balances)
    ch_energy_balances = xr.open_dataarray(path_to_ch_energy_balances)
    all_energy_balances = utils.merge_da([eurostat_energy_balances, ch_energy_balances], "annual_energy_balances")

    carrier_name_mapping = pd.read_csv(path_to_carrier_names, index_col=0)
    if sector == "commercial":
        # ASSUME: building heat demand from "other" sectors (incl. military, fisheries, agriculture)
        # is assigned to the commercial Euro-Calliope "sector"
        commercial_energy_balance = rename_carriers(
            slice_on_sector(all_energy_balances, "commercial"),
            carrier_name_mapping, "commercial"
        )
        other_energy_balance = rename_carriers(
            slice_on_sector(all_energy_balances, "other"),
            carrier_name_mapping, "other"
        ).reindex_like(commercial_energy_balance).fillna(0)

        sectoral_energy_balance = commercial_energy_balance + other_energy_balance

        assert (sectoral_energy_balance.groupby("carrier_name").sum(...) >= commercial_energy_balance.groupby("carrier_name").sum(...)).all()

    elif sector == "industry":
        sectoral_energy_balance = rename_carriers(
            eurostat_industry_subsectors_to_jrc_names(all_energy_balances, path_to_category_names),
            carrier_name_mapping, sector
        )
    else:
        sectoral_energy_balance = rename_carriers(
            slice_on_sector(all_energy_balances, sector),
            carrier_name_mapping, sector
        )

    sectoral_energy_balance = sectoral_energy_balance.sel(year=YEARS)

    sectoral_energy_balance.to_netcdf(path_to_output)


def eurostat_industry_subsectors_to_jrc_names(all_energy_balances, category_name_mapping):
    return utils.rename_and_groupby(
        all_energy_balances, category_name_mapping,
        dim_name="cat_code", new_dim_name="cat_name"
    )


def slice_on_sector(all_energy_balances, sector):
    return (
        all_energy_balances
        .sel(cat_code=EUROSTAT_SECTOR_MAPPING[sector])
        .sum("cat_code", min_count=1, keep_attrs=True)
    )


def rename_carriers(sectoral_energy_balance, carrier_name_mapping, sector):
    sectoral_energy_balance = utils.rename_and_groupby(
        sectoral_energy_balance, carrier_name_mapping[f"{sector}_carrier_name"].dropna().to_dict(),
        dim_name="carrier_code", new_dim_name="carrier_name"
    )

    return sectoral_energy_balance



if __name__ == "__main__":
    sectoral_energy_balances(
        path_to_eurostat_energy_balances=snakemake.input.eurostat_energy_balances,
        path_to_ch_energy_balances=snakemake.input.ch_energy_balances,
        sector=snakemake.wildcards.building_sector,
        path_to_carrier_names=snakemake.params.carrier_names,
        path_to_category_names=snakemake.params.category_names,
        path_to_output=snakemake.output[0]
    )
