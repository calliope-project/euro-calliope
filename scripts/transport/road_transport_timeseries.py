import pandas as pd
import datetime


def create_road_transport_demand_timeseries(
    road_distance_path, road_distance_timeseries_out_path
):
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
            expanded_data.append([
                current_date.strftime('%Y-%m-%d %H:%M:%S'),
                vehicle,
                country,
                hourly_distance * conversion_factor  # here I can just multiply this with a conversion factor
                #to get the desired energy demand
            ])

    expanded_df = pd.DataFrame(expanded_data, columns=["timestamp", "vehicle_type", "country_code", "energy_demand"])

    expanded_df.to_csv(road_distance_timeseries_out_path, index=False)


if __name__ == "__main__":
    create_road_transport_demand_timeseries(
        road_distance_path=snakemake.input.road_distance_path,
        road_distance_timeseries_out_path=snakemake.output.road_distance_timeseries_out_path,
    )
