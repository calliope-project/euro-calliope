import numpy as np
import pandas as pd
# import snakemake
from utils import jrc_idees_parser as jrc

CAT_NAME_STEEL = "Iron and steel"

H2_LHV_KTOE = 2.863  # 0.0333 TWh/kt LHV -> 2.863ktoe/kt
RECYCLED_STEEL = 0.5  # 50% H-DRI Iron, 50% scrap steel
H2_TO_STEEL = (1 - RECYCLED_STEEL) * 0.05  # 0.05t_h2/t_steel in H-DRI
HDRI_CONSUMPTION = 0.0116  # H-DRI: 135kWh_e/t = 0.0116ktoe/kt
MWH_PER_T_TO_KTOE_PER_KT = 0.08598  # 1 MWh/t -> 1 GWh/kt -> 0.08598 ktoe/kt
METHANOL_LHV_KTOE = 0.476  # 19.915 MJ/kg LHV -> 19915000MJ/kt -> 0.476ktoe/kt

def get_steel_iron_demand(path_jrc_energy, path_jrc_production):
    # Ensure we only operate with the relevant sector's data
    energy_df = pd.read_csv(path_jrc_energy, index_col=[0, 1, 2, 3, 4, 5, 6])
    energy_df = energy_df.xs(CAT_NAME_STEEL, level="cat_name", drop_level=False)
    prod_df = pd.read_csv(path_jrc_production, index_col=[0, 1, 2, 3])
    prod_df = prod_df.xs(CAT_NAME_STEEL, level="cat_name", drop_level=False)

    demand = energy_df.xs('demand').sum(level=['section', 'subsection', 'country_code', 'cat_name', 'unit'])

    # TODO: acquiring all the demands below should be done in a utils function

    # if it can be met by electricity (exclusively or otherwise),
    # then it's an end-use electricity demand
    electrical_consumption = (
        jrc.get_carrier_demand('Electricity', demand, energy_df)
        .assign(carrier='electricity').set_index('carrier', append=True)
    )
    # If it can only be met by natural gas (steam heating) then it's natural gas
    nat_gas_consumption = (
        jrc.get_carrier_demand('Natural gas (incl. biogas)', demand, energy_df)
        .drop(electrical_consumption.droplevel('carrier').index, errors='ignore')
        .assign(carrier='methane').set_index('carrier', append=True)
    )
    # If it can only be met by diesel (backup generators) then it's diesel
    diesel_consumption = (
        jrc.get_carrier_demand('Diesel oil (incl. biofuels)', demand, energy_df)
        .drop(nat_gas_consumption.droplevel('carrier').index, errors='ignore')
        .drop(electrical_consumption.droplevel('carrier').index, errors='ignore')
        .assign(carrier='diesel').set_index('carrier', append=True)
    )

    # TODO: space heat demand should be taken here, in a generic form

    steel_energy_consumption = get_steel_energy_consumption(energy_df, prod_df)
    pass


def get_steel_energy_consumption(energy_df, prod_df):
    """
    Calculates energy consumption in the iron and steel industry based on expected
    change in process to avoid fossil feedstocks. All process specific energy consumption
    (energy/t_steel) is based on the Electric Arc process (EAF), except sintering, which
    will be required for iron ore processed using H-DRI, but is not required by EAF.

    This function does the following:
    1. Finds all the specific consumption values by getting
        a. process energy demand / produced steel => specific demand
        b. process electrical demand / electrical consumption => electrical efficiency
        c. specific demand / electricial efficiency => specific electricity consumption
    2. Gets total process specific electricity consumption by adding specific consumptions
    for direct electric processes, EAF, H-DRI, smelting, sintering, refining, and finishing
    3. Gets specific hydrogen consumption for all countries that will process iron ore
    4. Gets specific space heat demand based on demand associated with EAF plants
    5. Gets total demand for electricity, hydrogen, and space heat by multiplying specific
    demand by total steel production (by both EAF and BF-BOF routes).
    """

    # sintering/pelletising
    sintering_specific_consumption = jrc.get_specific_electricity_consumption(
        'Integrated steelworks', 'Steel: Sinter/Pellet making',
        energy_df, prod_df
    )

    # smelters
    eaf_smelting_specific_consumption = jrc.get_specific_electricity_consumption(
        'Electric arc', 'Steel: Smelters',
        energy_df, prod_df
    )

    # EAF
    eaf_specific_consumption = jrc.get_specific_electricity_consumption(
        'Electric arc', 'Steel: Electric arc',
        energy_df, prod_df
    )

    # Rolling & refining
    refining_specific_consumption = jrc.get_specific_electricity_consumption(
        'Electric arc', 'Steel: Furnaces, Refining and Rolling',
        energy_df, prod_df
    )
    finishing_specific_consumption = jrc.get_specific_electricity_consumption(
        'Electric arc', 'Steel: Products finishing',
        energy_df, prod_df
    )

    # Auxiliaries (lighting, motors, etc.)
    auxiliary_specific_consumption = jrc.get_auxiliary_electricity_consumption(
        'Electric arc', energy_df, prod_df
    )

    # Total electricity consumption
    ## If the country produces steel from Iron ore (assuming 50% recycling):
    ### sintering/pelletising * iron_ore_% + smelting * recycled_steel_% + H-DRI + EAF + refining/rolling + finishing + auxiliaries
    ## If the country only recycles steel:
    ### smelting + EAF + refining/rolling + finishing + auxiliaries

    total_specific_consumption = (
        sintering_specific_consumption.mul(1 - RECYCLED_STEEL).add(
            eaf_smelting_specific_consumption
            # if no sintering, this country/year recycles 100% of steel
            .where(sintering_specific_consumption == 0)
            # if there is some sintering, we update smelting consumption to equal our assumed 2050 recycling rate and add weighted H-DRI consumption to process the remaining iron ore
            .fillna(eaf_smelting_specific_consumption.mul(RECYCLED_STEEL).add(HDRI_CONSUMPTION))
        )
        .add(eaf_specific_consumption)
        .add(refining_specific_consumption)
        .add(finishing_specific_consumption)
        .add(auxiliary_specific_consumption)
    )
    # In case our model now says a country does produce steel,
    # we give them the average of energy consumption of all other countries
    total_specific_consumption = (
        total_specific_consumption
        .where(total_specific_consumption > 0)
        .fillna(total_specific_consumption.mean())
        .assign(carrier='electricity')
        .set_index('carrier', append=True)
    )

    # Hydrogen consumption for H-DRI, only for those country/year combinations that handle iron ore
    # and don't recycle all their steel
    h2_specific_consumption = H2_LHV_KTOE * H2_TO_STEEL
    total_specific_h2_consumption = (
        total_specific_consumption
        .where(sintering_specific_consumption > 0)
        .fillna(0)
        .where(lambda x: x == 0)
        .fillna(h2_specific_consumption)
        .rename(index={'electricity': 'hydrogen'})
    )
    total_specific_consumption = total_specific_consumption.append(total_specific_h2_consumption)
    # Space heat
    space_heat_specific_demand = (
        energy_df
        .xs(('demand', 'Electric arc', 'Low enthalpy heat'))
        .div(prod_df.xs('Electric arc').droplevel('unit'))
        .assign(carrier='space_heat').set_index('carrier', append=True)
        .sum(level=total_specific_consumption.index.names)
        .rename(index={'ktoe': 'ktoe/kt'})
    )
    total_specific_consumption = total_specific_consumption.append(space_heat_specific_demand)

    steel_consumption = (
        total_specific_consumption
        .mul(prod_df.xs('Iron and steel', level='cat_name').sum(level='country_code'), level='country_code')
        .rename(index={'ktoe/kt': 'ktoe'})
    )

    return steel_consumption


if __name__ == "__main__":
    get_steel_iron_demand(
        path_jrc_energy=snakemake.input.path_jrc_energy,
        path_jrc_production=snakemake.input.path_jrc_production
    )


