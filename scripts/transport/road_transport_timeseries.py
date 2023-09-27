import pandas as pd
import datetime


def create_road_transport_demand_timeseries(
    road_distance_path, bau_electricity_path, first_year, final_year,
    power_scaling_factor, ldv_conversion_factor, hdv_conversion_factor, coaches_and_buses_conversion_factor,
    passenger_cars_conversion_factor, powered_2_wheelers_conversion_factor, light_duty_vehicles_timeseries_out_path,
    heavy_duty_vehicles_timeseries_out_path, coaches_and_buses_timeseries_out_path, passenger_cars_timeseries_out_path,
    powered_2_wheelers_timeseries_out_path, light_duty_vehicles_bau_timeseries_out_path,
    coaches_and_buses_bau_timeseries_out_path, passenger_cars_bau_timeseries_out_path
):
    # Read annual road transport distance into panda dataframe
    df = pd.read_csv(road_distance_path)

    # Read annual road transport bau electricity into panda dataframe
    df_bau = pd.read_csv(bau_electricity_path)

    # Assumption: Every hour of the day the same amount of distance gets driven.
    # Calculate the distance travelled per hour and add new column for it into road transport distance dataframe
    df["distance_per_hour"] = df["0"] / 8760

    # Assumption: Every hour of the day the same amount of distance gets driven.
    # Calculate the bau electricity per hour and add new column for it into road transport bau electricity dataframe
    df_bau["electricity_per_hour"] = df_bau["0"] / 8760

    # Process the road transport distance timeseries for each vehicle type
    process_timeseries("Light duty vehicles", power_scaling_factor, ldv_conversion_factor, first_year,
                       final_year, df, light_duty_vehicles_timeseries_out_path, False)
    process_timeseries("Heavy duty vehicles", power_scaling_factor, hdv_conversion_factor, first_year,
                       final_year, df, heavy_duty_vehicles_timeseries_out_path, False)
    process_timeseries("Motor coaches, buses and trolley buses", power_scaling_factor,
                       coaches_and_buses_conversion_factor, first_year, final_year, df,
                       coaches_and_buses_timeseries_out_path, False)
    process_timeseries("Passenger cars", power_scaling_factor, passenger_cars_conversion_factor, first_year,
                       final_year, df, passenger_cars_timeseries_out_path, False)
    process_timeseries("Powered 2-wheelers", power_scaling_factor, powered_2_wheelers_conversion_factor,
                       first_year, final_year, df, powered_2_wheelers_timeseries_out_path, False)

    # Process the road transport bau electricity timeseries for each vehicle type
    process_timeseries("Light duty vehicles", power_scaling_factor, 1, first_year, final_year, df_bau,
                       light_duty_vehicles_bau_timeseries_out_path, True)
    process_timeseries("Motor coaches, buses and trolley buses", power_scaling_factor, 1, first_year,
                       final_year, df_bau, coaches_and_buses_bau_timeseries_out_path, True)
    process_timeseries("Passenger cars", power_scaling_factor, 1, first_year, final_year, df_bau,
                       passenger_cars_bau_timeseries_out_path, True)


def process_timeseries(vehicle, power_scaling_factor, conversion_factor, first_year, final_year, df,
                       output_path, bau):
    # Filter data for the current vehicle type
    vehicle_df = df[df["vehicle_type"] == vehicle]

    # Create a df with index as timestamps and columns as countries.
    # The Dataframe is initially filled with zeros.
    timeseries_df = pd.DataFrame(index=[datetime.datetime(year, 1, 1) + datetime.timedelta(hours=hour) for year in
                                 range(first_year, final_year+1) for hour in range(8760)],
                                 columns=vehicle_df["country_code"].unique()).fillna(0)

    # Populate the timeseries with either distance per hour for road transport distance or electricity per hour for
    # road transport bau electricity
    for _, row in vehicle_df.iterrows():
        # For each row in vehicle_df, get the year and country code
        curr_year = row['year']
        country_code = row['country_code']

        # Create mask to find corresponding rows in timeseries_df
        mask = timeseries_df.index.year == curr_year

        # BAU case: set cells in timeseris_df with electricty per hour, multiply this value with the
        # power_scaling_factor out of the config, and a conversion factor which is 1, since we already have electricity.
        if bau:
            timeseries_df.loc[mask, country_code] = row["electricity_per_hour"] * conversion_factor * power_scaling_factor
        # Set cells in timeseris_df with distance per hour multiplied with a conversion factor depending on the vehicle
        # type to get from mio km to Mw(h), after that apply the power scaling factor to get to 100 GW(h) from Mw(h)
        # Finally make that value negative, since calliope works with negative values for demands.
        else:
            timeseries_df.loc[mask, country_code] = row["distance_per_hour"] * conversion_factor * power_scaling_factor * -1
    # Save the created timeseries to a csv file.
    timeseries_df.to_csv(output_path, index_label="utc-timestamp")


if __name__ == "__main__":
    create_road_transport_demand_timeseries(
        # input paths to total road distance and bau electricity
        road_distance_path=snakemake.input.road_distance_path,
        bau_electricity_path=snakemake.input.bau_electricity_path,
        # power scaling factor used to get from MW(h) to GW(h)
        power_scaling_factor=snakemake.params.power_scaling_factor,  # temporal scope
        first_year=snakemake.params.first_year, final_year=snakemake.params.final_year,
        # conversion factors for all vehicle types
        ldv_conversion_factor=snakemake.params.ldv_conversion_factor,
        hdv_conversion_factor=snakemake.params.hdv_conversion_factor,
        coaches_and_buses_conversion_factor=snakemake.params.coaches_and_buses_conversion_factor,
        passenger_cars_conversion_factor=snakemake.params.passenger_cars_conversion_factor,
        powered_2_wheelers_conversion_factor=snakemake.params.powered_2_wheelers_conversion_factor,
        # output paths for vehicle type electricity timeseries
        light_duty_vehicles_timeseries_out_path=snakemake.output.light_duty_vehicles_timeseries_out_path,
        heavy_duty_vehicles_timeseries_out_path=snakemake.output.heavy_duty_vehicles_timeseries_out_path,
        coaches_and_buses_timeseries_out_path=snakemake.output.coaches_and_buses_timeseries_out_path,
        passenger_cars_timeseries_out_path=snakemake.output.passenger_cars_timeseries_out_path,
        powered_2_wheelers_timeseries_out_path=snakemake.output.powered_2_wheelers_timeseries_out_path,
        # output paths for bau electricity timeseries per vehicle type
        light_duty_vehicles_bau_timeseries_out_path=snakemake.output.light_duty_vehicles_bau_timeseries_out_path,
        coaches_and_buses_bau_timeseries_out_path=snakemake.output.coaches_and_buses_bau_timeseries_out_path,
        passenger_cars_bau_timeseries_out_path=snakemake.output.passenger_cars_bau_timeseries_out_path, )
