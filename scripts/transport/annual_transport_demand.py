from enum import Enum

import pandas as pd

from eurocalliopelib import utils

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


class FinalConsumption(str, Enum):
    ROAD_TRANSPORT = "FC_TRA_ROAD_E"
    OTHER_SECTORS = "FC_OTH_NSP_E"
    AGRICULTURE_AND_FORESTRY = "FC_OTH_AF_E"


class Carrier(str, Enum):
    AVIATION_GASOLINE = "O4651"
    GASOLINE_JET_FUEL = "O4653"
    KEROSENE_JET_FUEL = "O4661XR5230B"
    OTHER_KEROSENE = "O4669"
    OIL_AND_PETROLEUM = "O4000XBIO"


AVIATION_CARRIERS = [
    Carrier.AVIATION_GASOLINE,
    Carrier.GASOLINE_JET_FUEL,
    Carrier.KEROSENE_JET_FUEL,
    Carrier.OTHER_KEROSENE
]


def get_all_distance_efficiency(energy_balance: pd.Series,
                                other_transport_road: pd.Series, road_energy: pd.Series, road_distance: pd.Series,
                                fill_missing_values: dict[str, str]
                                ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Add transport energy demand from agriculture and 'not elsewhere specified' (military) (OTHER_TRANSPORT_ROAD)
    transport_energy_balance = (
        energy_balance
        .xs(FinalConsumption.ROAD_TRANSPORT)
        .unstack('carrier_code')
        .groupby(CARRIERS, axis=1).sum(min_count=1)
        .rename_axis(columns='carrier')
        .droplevel('unit')
        .stack()
        .add(other_transport_road, fill_value=0)
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
        road_energy.div(road_energy.groupby(level=['carrier', 'country_code', 'year']).sum()),
        fill_missing_values
    )

    # Energy consumption per transport mode by mapping transport mode
    # carrier contributions to total carrier consumption
    transport_energy_per_mode = carrier_contribution.mul(transport_energy_balance).dropna()

    # Distance per unit energy consumed per transport mode according to JRC IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    transport_efficiency = fill_missing_countries_and_years(
        road_distance
        .where(road_distance > 0)
        .div(
            road_energy
            .where(road_energy > 0)
            .groupby(level=['country_code', 'vehicle_subtype', 'section', 'vehicle_type', 'year'])
            .sum()
        ),
        fill_missing_values
    )

    # Distance travelled per transport mode, including years 2016-2018,
    # based on JRC IDEES transport efficiency (2015 data for 2016-2018)
    transport_distance_all_years = (
        transport_energy_per_mode
        .groupby(level=['country_code', 'vehicle_subtype', 'section', 'vehicle_type', 'year'])
        .sum()
        .mul(transport_efficiency)
    )

    # Use the distance traveled based on Eurostat to fill in blanks in JRC data
    aligned_dfs = (
        transport_distance_all_years
        .reorder_levels(road_distance.index.names)
        .align(road_distance)
    )
    total_transport_distance = (
        aligned_dfs[1]
        .fillna(aligned_dfs[0])
        .groupby(level=['vehicle_type', 'vehicle_subtype', 'country_code', 'year'])
        .sum()
    )

    return total_transport_distance, transport_efficiency, transport_energy_per_mode


def get_all_vehicles(total_road_distance: pd.DataFrame, road_vehicles: pd.Series,
                     fill_missing_values: dict[str, str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    total_vehicle_distance = fill_missing_countries_and_years(
        road_vehicles.div(road_vehicles),
        fill_missing_values
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


def fill_missing_countries_and_years(jrc_data: pd.DataFrame, fill_missing_values: dict[str, str]) -> pd.DataFrame:
    jrc_data = jrc_data.unstack('country_code')
    try:
        jrc_data = jrc_data.loc[:, "value"]
    except KeyError:
        pass # it's fine. Just checking there is no MultiIndex in the columns
    for country, neighbors in fill_missing_values.items():
        jrc_data = jrc_data.assign(**{country: jrc_data[neighbors].mean(axis=1)})

    jrc_data = jrc_data.stack().unstack('year')
    jrc_data = jrc_data.assign(
        **{str(year): jrc_data[2015] for year in range(2016, 2019)}
    )
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()


def road_efficiency_cleanup(road_efficiency: pd.DataFrame, efficiency_quantile: float) -> pd.DataFrame:
    return (
        (1 / road_efficiency.xs(2015, level='year'))
        .unstack(['vehicle_type', 'vehicle_subtype'])
        .quantile(efficiency_quantile)
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


def get_historic_road_electricity_consumption(road_historic_consumption: pd.DataFrame) -> pd.DataFrame:
    return (
        road_historic_consumption
        .groupby(level=['carrier', 'vehicle_type', 'country_code', 'year'])
        .sum()
        .xs('electricity')
    )


if __name__ == "__main__":
    energy_balances = pd.read_csv(
        snakemake.input.energy_balances,
        index_col=["cat_code", "carrier_code", "unit", "country", "year"]
    ).squeeze()
    road_energy = pd.read_csv(
        snakemake.input.jrc_road_energy,
        index_col=["section", "vehicle_type", "vehicle_subtype", "carrier", "country_code", "year"]
    ).squeeze()
    road_distance = pd.read_csv(
        snakemake.input.jrc_road_distance,
        index_col=["section", "vehicle_type", "vehicle_subtype", "country_code", "year"]
    ).squeeze()
    road_vehicles = pd.read_csv(
        snakemake.input.jrc_road_vehicles,
        index_col=["section", "vehicle_type", "vehicle_subtype", "country_code", "year"]
    )
    # Used to add transport energy demand from agriculture and 'not elsewhere specified' (military)
    # ASSUME: agriculture oil use goes to 'road' transport demand;
    # 'not elsewhere specified' oil use goes predominantly to 'road' transport, except kerosene which goes to aviation
    other_transportation_aviation = (  # all kerosene from the military destined for aviation
        energy_balances
        .loc[idx[FinalConsumption.OTHER_SECTORS, AVIATION_CARRIERS, :, :, :]]
        .groupby(level=["country", "year"])
        .sum()
    )
    other_transport_road = (  # i.e. all oil that isn't destined for aviation
        energy_balances
        .loc[idx[[FinalConsumption.AGRICULTURE_AND_FORESTRY, FinalConsumption.OTHER_SECTORS], Carrier.OIL_AND_PETROLEUM, :, :, :]]
        .groupby(level=["country", "year"])
        .sum()
        .sub(other_transportation_aviation, fill_value=0)  # remove fuel use assumed for aviation
        .to_frame('diesel')
        .rename_axis(columns='carrier')
        .stack()
    )

    fill_missing_values = snakemake.params.fill_missing_values

    # Calculate total road distance, road efficiency and historically electrified road consumption
    total_road_distance, road_efficiency, road_historically_electrified_consumption = get_all_distance_efficiency(
        energy_balance=energy_balances,
        other_transport_road=other_transport_road,
        road_energy=road_energy,
        road_distance=road_distance,
        fill_missing_values=fill_missing_values
    )

    # Calculate total road vehicles and aggregate total road distance
    total_road_vehicles, total_road_distance = get_all_vehicles(
        total_road_distance,
        road_vehicles,
        fill_missing_values
    )

    # Some cleanup that's specific to road data for road efficiency
    road_efficiency = road_efficiency_cleanup(road_efficiency, snakemake.params.efficiency_quantile)

    # Extract historical electricity consumption
    road_electricity_historic = get_historic_road_electricity_consumption(road_historically_electrified_consumption)

    # Create CSV Files for calculated data
    total_road_distance.rename("value").to_csv(snakemake.output.distance)
    total_road_vehicles.rename("value").to_csv(snakemake.output.vehicles)
    road_efficiency.rename("value").to_csv(snakemake.output.efficiency)
    road_electricity_historic.rename("value").to_csv(snakemake.output.road_electricity_historic)
