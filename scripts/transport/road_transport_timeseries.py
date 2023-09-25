import pandas as pd
import datetime


def create_road_transport_demand_timeseries(road_distance_path, road_distance_timeseries_out_path):
    conversion_factor = 1
    road_distance_df = pd.read_csv(road_distance_path)

    expanded_data = []

    for _, row in road_distance_df.iterrows():
        vehicle = row["vehicle_type"]
        country = row["country_code"]
        print(country)
        year = row["year"]
        yearly_distance = row["0"]
        start_date = datetime.datetime(year, 1, 1)
        # calculate the distance traveled per hour
        # assumption each hour a day the same amount of kilometers get driven
        hourly_distance = yearly_distance / 8760
        for hour in range(8760):
            current_date = start_date + datetime.timedelta(hours=hour)
            expanded_data.append(
                [current_date.strftime('%Y-%m-%d %H:%M:%S'), vehicle, country, hourly_distance * conversion_factor
                 # here I can just multiply this with a conversion factor
                 # to get the desired energy demand
                 # TODO: bau muss positiv sein und andere negativ!
                 ])

    expanded_df = pd.DataFrame(expanded_data, columns=["timestamp", "vehicle_type", "country_code", "energy_demand"])

    expanded_df.to_csv(road_distance_timeseries_out_path, index=False)


def process_timeseries_for_vehicle(vehicle, power_scaling_factor, conversion_factor, first_year, final_year, df,
                                   output_path, bau):
    print("in function")
    # Filter data for the current vehicle type
    vehicle_df = df[df["vehicle_type"] == vehicle]
    print(vehicle_df)

    # Create a df with index as timestamps and columns as countries.
    # The table is initially filled with zeros.
    timeseries_df = pd.DataFrame(index=[datetime.datetime(year, 1, 1) + datetime.timedelta(hours=hour) for year in
                                      range(first_year, final_year + 1) for hour in range(8760)],
        columns=vehicle_df["country_code"].unique()).fillna(0)
    print("this is the conversion_factor: ")
    print(conversion_factor)
    print("this is the power_scaling_factor: ")
    print(power_scaling_factor)
    # Populate the timeseries
    for _, row in vehicle_df.iterrows():
        # For each row in vehicle_df, get the year and country code
        curr_year = row['year']
        print("curr_year = ")
        print(curr_year)
        country_code = row['country_code']
        print("country code = ")
        print(country_code)

        # Find the corresponding rows in the pivot_table based on the year
        mask = timeseries_df.index.year == curr_year

        # Update the corresponding cell in the pivot_table
        if bau:
            timeseries_df.loc[mask, country_code] = row["electricity_per_hour"] * conversion_factor * power_scaling_factor
        else:
            timeseries_df.loc[mask, country_code] = row[
                                                      "distance_per_hour"] * conversion_factor * power_scaling_factor * -1
    timeseries_df.to_csv(output_path, index_label="utc-timestamp")
    print(vehicle + " timeseries is done")


def create_road_transport_demand_timeseries_v2(road_distance_path, bau_electricity_path, first_year, final_year,
    power_scaling_factor, ldv_conversion_factor, hdv_conversion_factor, coaches_and_buses_conversion_factor,
    passenger_cars_conversion_factor, powered_2_wheelers_conversion_factor, light_duty_vehicles_timeseries_out_path,
    heavy_duty_vehicles_timeseries_out_path, coaches_and_buses_timeseries_out_path, passenger_cars_timeseries_out_path,
    powered_2_wheelers_timeseries_out_path, light_duty_vehicles_bau_timeseries_out_path,
    coaches_and_buses_bau_timeseries_out_path, passenger_cars_bau_timeseries_out_path):

    df = pd.read_csv(road_distance_path)

    df_bau = pd.read_csv(bau_electricity_path)
    print(df_bau)

    # Calculate the distance travelled per hour
    df["distance_per_hour"] = df["0"] / 8760

    # Calculate the bau electricity per hour
    df_bau["electricity_per_hour"] = df_bau["0"] / 8760
    print(df_bau)
    # Process the timeseries for each vehicle type
    process_timeseries_for_vehicle("Light duty vehicles", power_scaling_factor, ldv_conversion_factor, first_year,
                                   final_year, df, light_duty_vehicles_timeseries_out_path, False)
    process_timeseries_for_vehicle("Heavy duty vehicles", power_scaling_factor, hdv_conversion_factor, first_year,
                                   final_year, df, heavy_duty_vehicles_timeseries_out_path, False)
    process_timeseries_for_vehicle("Motor coaches, buses and trolley buses", power_scaling_factor,
                                   coaches_and_buses_conversion_factor, first_year, final_year, df,
                                   coaches_and_buses_timeseries_out_path, False)
    process_timeseries_for_vehicle("Passenger cars", power_scaling_factor, passenger_cars_conversion_factor, first_year,
                                   final_year, df, passenger_cars_timeseries_out_path, False)
    process_timeseries_for_vehicle("Powered 2-wheelers", power_scaling_factor, powered_2_wheelers_conversion_factor,
                                   first_year, final_year, df, powered_2_wheelers_timeseries_out_path, False)

    # process timeseries for each vehicle type for bau
    process_timeseries_for_vehicle("Light duty vehicles", power_scaling_factor, 1, first_year, final_year, df_bau,
                                   light_duty_vehicles_bau_timeseries_out_path, True)
    process_timeseries_for_vehicle("Motor coaches, buses and trolley buses", power_scaling_factor, 1, first_year,
                                   final_year, df_bau, coaches_and_buses_bau_timeseries_out_path, True)
    process_timeseries_for_vehicle("Passenger cars", power_scaling_factor, 1, first_year, final_year, df_bau,
                                   passenger_cars_bau_timeseries_out_path, True)


if __name__ == "__main__":
    create_road_transport_demand_timeseries_v2(  # input paths to total road distance and bau electricity
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
