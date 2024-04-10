import eurocalliopelib.utils as ec_utils
import pandas as pd
from utils import formatting
from utils import jrc_idees_parser as jrc


H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
MWH_PER_T_TO_KTOE_PER_KT = 0.08598  # 1 MWh/t -> 1 GWh/kt -> 0.08598 ktoe/kt
METHANOL_LHV_KTOE = 0.476  # 19.915 MJ/kg LHV -> 19915000MJ/kt -> 0.476ktoe/kt

CAT_NAME_CHEM = "Chemicals Industry"


def get_chem_demand_df(
    year_range: list,
    path_energy_balances: str,
    path_cat_names: str,
    path_carrier_names: str,
    path_jrc_energy: str,
    path_jrc_production: str,
    path_output: str = "",
):
    # -------------------------------------------------------------------------
    # Prepare data files
    # -------------------------------------------------------------------------
    energy_df = pd.read_csv(path_jrc_energy, index_col=[0, 1, 2, 3, 4, 5, 6])
    prod_df = pd.read_csv(path_jrc_production, index_col=[0, 1, 2, 3])
    cat_names_df = pd.read_csv(path_cat_names, header=0, index_col=0)
    carrier_names_df = pd.read_csv(path_carrier_names, header=0, index_col=0)
    # Ensure dataframes only have data specific to this industry
    cat_names_df = cat_names_df[cat_names_df["jrc_idees"] == CAT_NAME_CHEM]
    energy_df = energy_df.xs(CAT_NAME_CHEM, level="cat_name", drop_level=False)
    prod_df = prod_df.xs(CAT_NAME_CHEM, level="cat_name", drop_level=False)
    demand = energy_df.xs('demand').sum(level=['section', 'subsection', 'country_code', 'cat_name', 'unit'])
    energy_balances_df = pd.read_csv(
        path_energy_balances, index_col=[0, 1, 2, 3, 4], squeeze=True
    )
    # if it can be met by electricity (exclusively or otherwise),
    # then it's an end-use electricity demand
    electrical_consumption = (
        jrc.get_carrier_demand('Electricity', demand, energy_df)
        .assign(carrier='electricity').set_index('carrier', append=True)
    )
    # -------------------------------------------------------------------------
    # Process data
    # -------------------------------------------------------------------------
    chemicals_energy_consumption = process_chemicals_energy_consumption(electrical_consumption, prod_df, demand)
    # -------------------------------------------------------------------------
    # Format the final output
    # -------------------------------------------------------------------------
    chemicals_energy_consumption.columns = chemicals_energy_consumption.columns.astype(
        int
    ).rename("year")
    filled_consumption_df = formatting.fill_missing_data(
        energy_balances_df,
        cat_names_df,
        carrier_names_df,
        chemicals_energy_consumption,
        year_range,
    )

    units = filled_consumption_df.index.get_level_values("unit")
    filled_consumption_df.loc[units == "ktoe"] = filled_consumption_df.loc[
        units == "ktoe"
    ].apply(ec_utils.ktoe_to_twh)
    filled_consumption_df = filled_consumption_df.rename({"ktoe": "twh"}, level="unit")
    filled_consumption_df.index = filled_consumption_df.index.set_names(
        "subsector", level="cat_name"
    )
    filled_consumption_df = filled_consumption_df.stack()

    if path_output:
        filled_consumption_df.reorder_levels(formatting.LEVEL_ORDER).to_csv(path_output)

    return filled_consumption_df



def process_chemicals_energy_consumption(electrical_consumption, prod_df, demand):
    """
    We remove feedstock and steam-based processing from "basic chemicals" production,
    which will be replaced by direct H2 provision.
    All other demand (non-basic chemicals & basic chemicals electricity) is directly passed over.
    """
    chem_electricity_consumption = (
        electrical_consumption
        .drop('Low enthalpy heat', level='subsection')
        .xs('Chemicals Industry', level='cat_name')
        .sum(level=['country_code', 'unit', 'carrier'])
    )

    h2_demand = {  # t/t
        'Ethylene': 0,
        'Propylene': 0,
        'BTX': 0,
        'Ammonia': 0.178,
        'Methanol': 0.189,
    }
    methanol_demand = {
        'Ethylene': 2.83,
        'Propylene': 2.83,
        'BTX': 4.3,
        'Ammonia': 0,
        'Methanol': 1
    }
    co2_demand = {  # tCO2/t
        'Ethylene': 0,
        'Propylene': 0,
        'BTX': 0,
        'Ammonia': 0.112,  # inc. 0.35t_urea/t_ammonia
        'Methanol': 1.373,
    }

    energy_demand = {  # MWh/t
        'Ethylene': 1.4,  # MTO process
        'Propylene': 1.4,  # MTO process
        'BTX': 1.4,  # MTA process
        'Ammonia': 2.05,  # inc. 0.35t_urea/t_ammonia
        'Methanol': 1.5,  # auxiliary demand, could just be assumed as already included in existing electricity demand
    }

    mass = {  # Bazzanella and Ausfelder, 2017
        'Ethylene': 21.7,
        'Propylene': 17,
        'BTX': 15.7,
        'Ammonia': 17,
        'Methanol': 2,
    }
    molar_mass = {
        'Ethylene': 28.05,
        'Propylene': 42.08,
        'BTX': 93,
        'Ammonia': 17.01,
        'Methanol': 32.04,
    }
    moles = {k: mass[k] / molar_mass[k] for k in mass.keys()}
    molar_ratio = {k: v / sum(moles.values()) for k, v in moles.items()}

    basic_chemicals_mass = prod_df.xs('Basic chemicals (kt ethylene eq.)').sum(level=['country_code'])
    basic_chemicals_moles = basic_chemicals_mass / molar_mass['Ethylene']
    masses = {
        k: basic_chemicals_moles * molar_ratio[k] * molar_mass[k]
        for k in molar_ratio.keys()
    }

    chem_h2_consumption = sum(masses[k] * h2_demand[k] * H2_LHV_KTOE for k in h2_demand.keys())
    chem_methanol_consumption = sum(masses[k] * methanol_demand[k] * METHANOL_LHV_KTOE for k in methanol_demand.keys())
    chem_co2_consumption = sum(masses[k] * co2_demand[k] for k in co2_demand.keys())
    chem_energy_consumption = sum(masses[k] * energy_demand[k] * MWH_PER_T_TO_KTOE_PER_KT for k in energy_demand.keys())
    # Space heat
    space_heat_demand = (
        demand
        .xs('Low enthalpy heat', level='subsection')
        .xs(CAT_NAME_CHEM, level='cat_name')
        .sum(level=['country_code', 'unit'])
    )

    chem_consumption = pd.concat(
        [chem_electricity_consumption.add(chem_energy_consumption).droplevel('carrier'),
         chem_methanol_consumption.assign(unit='ktoe').set_index('unit', append=True),
         chem_h2_consumption.assign(unit='ktoe').set_index('unit', append=True),
         chem_co2_consumption.mul(1e3).assign(unit='t').set_index('unit', append=True),
         space_heat_demand],

        names=['carrier'], keys=['electricity', 'methanol', 'hydrogen', 'co2', 'space_heat']
    )

    chem_consumption = chem_consumption.assign(cat_name='Chemicals Industry').set_index('cat_name', append=True)

    return chem_consumption



if __name__ == "__main__":
    get_chem_demand_df(
    year_range=snakemake.params.year_range,
    path_energy_balances=snakemake.input.path_energy_balances,
    path_cat_names=snakemake.input.path_cat_names,
    path_carrier_names=snakemake.input.path_carrier_names,
    path_jrc_energy=snakemake.input.path_jrc_energy,
    path_jrc_production=snakemake.input.path_jrc_production,
    path_output=snakemake.output.path_output,
    )