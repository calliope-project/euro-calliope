import xarray as xr

JRC_IDEES_CARRIER_MAPPING = {
    'Diesel oil engine': 'diesel',
    'Battery electric vehicles': 'electricity',
    'Domestic': 'diesel',
    'International': 'diesel',
    'Gasoline engine': 'petrol',
    'Plug-in hybrid electric': 'petrol'
}


def get_road_transport_demand(
    path_to_energy_balances, path_to_jrc_transport_energy_consumption, path_to_jrc_transport_distance_travelled,
    path_to_transport_distance_output, path_to_transport_current_electricity_consumption_output,
):
    energy_balances = xr.open_dataarray(path_to_energy_balances)
    energy_consumption = xr.open_dataarray(path_to_jrc_transport_energy_consumption)
    distance_travelled = xr.open_dataarray(path_to_jrc_transport_distance_travelled)

    current_consumption = get_transport_energy_per_mode(energy_balances, energy_consumption)
    current_transport_distance = get_transport_distance_per_mode(
        energy_consumption, distance_travelled, current_consumption
    )

    current_electricity_consumption = current_consumption.sum("sector").sel(carrier="electricity")

    current_transport_distance.to_csv(path_to_transport_distance_output)
    current_electricity_consumption.to_csv(path_to_transport_current_electricity_consumption_output)


def get_transport_energy_per_mode(energy_balances, energy_consumption):
    sum_dims = ["vehicle_type", "section"]
    if "vehicle_subtype" in energy_consumption.dims:
        sum_dims.append("vehicle_subtype")

    # contribution of each transport mode to carrier consumption from JRC_IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    carrier_contribution = fill_missing_countries_and_years(
        energy_consumption / energy_consumption.sum(sum_dims)
    )
    # Energy consumption per transport mode by mapping transport mode
    # carrier contributions to total carrier consumption
    return energy_balances * carrier_contribution


def get_transport_distance_per_mode(
    energy_consumption, distance_travelled, transport_energy_per_mode
):
    sum_dims = []
    if "vehicle_subtype" in energy_consumption.dims:
        sum_dims.append("vehicle_subtype")
    # Distance per unit energy consumed per transport mode according to JRC IDEES
    # 2016-2018 from 2015 data; non-JRC countries, based on neighbour data
    transport_efficiency = fill_missing_countries_and_years(
        distance_travelled.where(distance_travelled > 0) /
        energy_consumption.where(energy_consumption > 0).sum(sum_dims)
    )
    # Distance travelled per transport mode, including years 2016-2018,
    # based on JRC IDEES transport efficiency (2015 data for 2016-2018)
    transport_distance_all_years = (
        transport_energy_per_mode.sum(sum_dims) * transport_efficiency
    )
    # Use the distance traveled based on Eurostat to fill in blanks in JRC data

    total_transport_distance = (
        distance_travelled
        .fillna(transport_distance_all_years)
        .sum(sum_dims + ["sector"])
    )

    return total_transport_distance


def fill_missing_countries_and_years(jrc_data):
    """
    ASSUME:
    1. For the countries not covered by JRC-IDEES, use average of all
    neighbouring / nearby countries.
    2. For all years 2016-2018, use 2015 data (the most recent in JRC-IDEES)
    """
    jrc_data = jrc_data.unstack('country_code')
    balkan_countries = jrc_data[['BG', 'HR', 'HU', 'RO', 'EL']].mean(axis=1)
    nordic_countries = jrc_data[['SE', 'DK']].mean(axis=1)
    ch_neighbours = jrc_data[['DE', 'AT', 'FR', 'IT']].mean(axis=1)
    jrc_data = jrc_data.assign(
        AL=balkan_countries,
        BA=balkan_countries,
        ME=balkan_countries,
        MK=balkan_countries,
        RS=balkan_countries,
        NO=nordic_countries,
        IS=nordic_countries,
        CH=ch_neighbours,
    ).stack().unstack('year')
    # Cannot 'assign' with a numeric key, so have to stringify and then convert back to integer
    jrc_data = jrc_data.assign(
        **{str(i): jrc_data[2015] for i in range(2016, 2019)}
    )
    jrc_data.columns = jrc_data.columns.astype(int)
    return jrc_data.stack()


if __name__ == "__main__":
    get_road_transport_demand(
        energy_balances_path=snakemake.input.energy_balances,
        jrc_road_energy_path=snakemake.input.jrc_road_energy,
        jrc_road_distance_path=snakemake.input.jrc_road_distance,
        jrc_road_vehicles_path=snakemake.input.jrc_road_vehicles,
        vehicle_efficiency_percentile=snakemake.params.vehicle_efficiency_percentile,
        road_distance_out_path=snakemake.output.distance,
        road_vehicles_out_path=snakemake.output.vehicles,
        road_efficiency_out_path=snakemake.output.efficiency,
        road_bau_electricity_out_path=snakemake.output.road_bau_electricity,
    )
