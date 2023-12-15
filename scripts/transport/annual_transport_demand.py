import pandas as pd

from eurocalliopelib import utils
from typing import Tuple

CARRIERS = {
    'O4652XR5210B': 'petrol',
    'R5210B': 'biofuels',
    'O4671XR5220B': 'diesel',
    'R5220B': 'biofuels',
    'O4630': 'lpg',
    'G3000': 'natural_gas',
    'E7000': 'electricity'
}
idx = pd.IndexSlice

YEAR_RANGE = slice(2000, 2018)

# DATA: Eurostat Energy Balances, JRC IDEES processed road energy/distance/vehicles
ENERGY_BALANCES = utils.read_tdf(snakemake.input.energy_balances)
ROAD_ENERGY_DF = utils.read_tdf(snakemake.input.jrc_road_energy)
ROAD_DISTANCE_DF = utils.read_tdf(snakemake.input.jrc_road_distance)
ROAD_VEHICLES_DF = utils.read_tdf(snakemake.input.jrc_road_vehicles)

# Used to add transport energy demand from agriculture and 'not elsewhere specified' (military)
# ASSUME: agriculture oil use goes to 'road' transport demand;
# 'not elsewhere specified' oil use goes predominantly to 'road' transport, except kerosene which goes to aviation
OTHER_TRANSPORT_AVIATION = (  # all kerosene from the military destined for aviation
    ENERGY_BALANCES
    .loc[idx['FC_OTH_NSP_E', ['O4651', 'O4653', 'O4661XR5230B', 'O4669'], :, :, :]]
    .sum(level=['country', 'year'])
)
OTHER_TRANSPORT_ROAD = (  # i.e. all oil that isn't destined for aviation
    ENERGY_BALANCES
    .loc[idx[['FC_OTH_AF_E', 'FC_OTH_NSP_E'], 'O4000XBIO', :, :, :]]
    .sum(level=['country', 'year'])
    .sub(OTHER_TRANSPORT_AVIATION, fill_value=0)  # remove fuel use assumed for aviation
    .to_frame('diesel')
    .rename_axis(columns='carrier')
    .stack()
)

# All countries not contained in the JRC IDEES, but in spatial scope
FILL_MISSING_VALUES = snakemake.params.fill_missing_values
# Vehicle Efficiency
EFFICIENCY_QUANTILE = snakemake.params.efficiency_quantile


def get_transport_demand(
    road_distance_out_path: str,
    road_vehicles_out_path: str,
    road_efficiency_out_path: str,
    road_bau_electricity_out_path: str
):
    # Calculate total road distance, road efficiency and road BAU consumption
    total_road_distance, road_efficiency, road_bau_consumption = get_all_distance_efficiency(
        'FC_TRA_ROAD_E',
        'vehicle_subtype'
    )

    # Calculate total road vehicles and aggregate total road distance
    total_road_vehicles, total_road_distance = get_all_vehicles(total_road_distance)

    # Some cleanup that's specific to road data for road efficiency
    road_efficiency = road_efficiency_cleanup(road_efficiency)

    # Extract electricity BAU consumption
    road_electricity_bau = get_road_electricity_bau_consumption(road_bau_consumption)

    # Create CSV Files for calculated data
    total_road_distance.to_csv(road_distance_out_path)
    total_road_vehicles.to_csv(road_vehicles_out_path)
    road_efficiency.to_csv(road_efficiency_out_path)
    road_electricity_bau.to_csv(road_bau_electricity_out_path)


def get_all_distance_efficiency(cat_name: str, unique_dim: str) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Add transport energy demand from agriculture and 'not elsewhere specified' (military) (OTHER_TRANSPORT_ROAD)
    transport_energy_balance = (
        ENERGY_BALANCES
        .xs(cat_name)
        .unstack('carrier_code')
        .groupby(CARRIERS, axis=1).sum(min_count=1)
        .rename_axis(columns='carrier')
        .droplevel('unit')
        .stack()
        .add(OTHER_TRANSPORT_ROAD, fill_value=0)
        .apply(utils.tj_to_twh)
        .rename_axis(index=['country_code', 'year', 'carrier'])
        .unstack("year")
        .loc[:, YEAR_RANGE]
        .interpolate(axis=1, limit_direction="both")
        .stack()
    )

    # contribution of each transport mode to carrier consumption from JRC_IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    carrier_contribution = fill_missing_countries_and_years(
        ROAD_ENERGY_DF
        .div(ROAD_ENERGY_DF.sum(level=['carrier', 'country_code', 'year'])),
    )

    # Energy consumption per transport mode by mapping transport mode
    # carrier contributions to total carrier consumption
    transport_energy_per_mode = carrier_contribution.mul(transport_energy_balance).dropna()

    # Distance per unit energy consumed per transport mode according to JRC IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    transport_efficiency = fill_missing_countries_and_years(
        ROAD_DISTANCE_DF
        .where(ROAD_DISTANCE_DF > 0)
        .div(
            ROAD_ENERGY_DF
            .where(ROAD_ENERGY_DF > 0)
            .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        )
    )

    # Distance travelled per transport mode, including years 2016-2018,
    # based on JRC IDEES transport efficiency (2015 data for 2016-2018)
    transport_distance_all_years = (
        transport_energy_per_mode
        .sum(level=['country_code', unique_dim, 'section', 'vehicle_type', 'year'])
        .mul(transport_efficiency)
    )

    # Use the distance traveled based on Eurostat to fill in blanks in JRC data
    aligned_dfs = (
        transport_distance_all_years
        .reorder_levels(ROAD_DISTANCE_DF.index.names)
        .align(ROAD_DISTANCE_DF)
    )
    total_transport_distance = (
        aligned_dfs[1]
        .fillna(aligned_dfs[0])
        .sum(level=['vehicle_type', unique_dim, 'country_code', 'year'])
    )

    return total_transport_distance, transport_efficiency, transport_energy_per_mode


def get_all_vehicles(total_road_distance: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    total_vehicle_distance = fill_missing_countries_and_years(
        ROAD_VEHICLES_DF
        .div(ROAD_DISTANCE_DF)
    )
    total_road_vehicles = total_road_distance.mul(
        total_vehicle_distance
        .mean(level=['vehicle_type', 'vehicle_subtype', 'country_code', 'year'])
    )
    total_road_distance = total_road_distance.sum(
        level=['vehicle_type', 'country_code', 'year']
    )
    total_road_vehicles = total_road_vehicles.sum(
        level=['vehicle_type', 'country_code', 'year']
    )

    return total_road_vehicles, total_road_distance


def fill_missing_countries_and_years(jrc_data: pd.DataFrame) -> pd.DataFrame:
    jrc_data = jrc_data.unstack('country_code')
    for country, neighbors in FILL_MISSING_VALUES.items():
        jrc_data = jrc_data.assign(**{country: jrc_data[neighbors].mean(axis=1)})

    jrc_data = jrc_data.stack().unstack('year')
    jrc_data = jrc_data.assign(
        **{str(i): jrc_data[2015] for i in range(2016, 2019)}
    )
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()


def road_efficiency_cleanup(road_efficiency: pd.DataFrame) -> pd.DataFrame:
    return (
        # 25th percentile of 2015 efficiency in TWh/mio km
        (1 / road_efficiency.xs(2015, level='year'))
        .unstack(['vehicle_type', 'vehicle_subtype'])
        .quantile(EFFICIENCY_QUANTILE)
        .unstack(0)
        .groupby({
            'Diesel oil engine': 'diesel',
            'Battery electric vehicles': 'electricity',
            'Domestic': 'diesel',
            'International': 'diesel',
            'Gasoline engine': 'petrol',
            'Plug-in hybrid electric': 'petrol'
        }).mean()
        .stack()
    )


def get_road_electricity_bau_consumption(road_bau_consumption: pd.DataFrame) -> pd.DataFrame:
    return (
        road_bau_consumption
        .sum(level=['carrier', 'vehicle_type', 'country_code', 'year'])
        .xs('electricity')
    )


if __name__ == "__main__":
    get_transport_demand(
        road_distance_out_path=snakemake.output.distance,
        road_vehicles_out_path=snakemake.output.vehicles,
        road_efficiency_out_path=snakemake.output.efficiency,
        road_bau_electricity_out_path=snakemake.output.road_bau_electricity,
    )
