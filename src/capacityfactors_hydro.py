import pandas as pd
import geopandas as gpd
import xarray as xr
from shapely.geometry import Point

WGS_84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


def main(path_to_locations, path_to_capacities, path_to_hydro_stations, threshold, path_to_ror, path_to_reservoir):
    locations = gpd.read_file(path_to_locations).to_crs(WGS_84).set_index("id")
    capacities = pd.read_csv(path_to_capacities, index_col=0)
    plants = xr.open_dataset(path_to_hydro_stations)
    hror = time_series(
        plants=plants.sel(id=(plants.type == "HROR")),
        locations=locations,
        capacities=capacities["installed_capacity_hror_MW"]
    )
    hror.where(hror >= threshold, 0).to_csv(path_to_ror, header=True, index=True)
    hdam = time_series(
        plants=plants.sel(id=(plants.type == "HDAM")),
        locations=locations,
        capacities=capacities["installed_capacity_hdam_MW"]
    )
    hdam.where(hdam >= threshold, 0).to_csv(path_to_reservoir, header=True, index=True)


def time_series(plants, locations, capacities):
    plant_centroids = gpd.GeoDataFrame(
        crs=WGS_84,
        geometry=list(map(Point, zip(plants.lon, plants.lat))),
        index=plants.id
    )
    location_of_plant = gpd.sjoin(plant_centroids, locations, how="left", op='intersects')["index_right"]
    location_of_plant = location_of_plant[~location_of_plant.index.duplicated()].rename_axis(index="id")
    return (plants.inflow_MWh.groupby(location_of_plant.to_xarray())
                             .sum(dim="id")
                             .to_dataframe()
                             .reset_index()
                             .pivot(index="time", columns="index_right", values="inflow_MWh")
                             .reindex(columns=locations.index, fill_value=0)
                             .div(capacities)
                             .rename(columns=lambda col: col.replace(".", "-")))


if __name__ == "__main__":
    main(
        path_to_locations=snakemake.input.locations,
        path_to_capacities=snakemake.input.capacities,
        path_to_hydro_stations=snakemake.input.stations,
        threshold=snakemake.params.threshold,
        path_to_ror=snakemake.output.ror,
        path_to_reservoir=snakemake.output.reservoir
    )
