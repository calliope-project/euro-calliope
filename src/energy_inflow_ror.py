import geopandas as gpd
import xarray as xr
from shapely.geometry import Point

WGS_84 = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"


def energy_inflow_per_location(path_to_locations, path_to_hydro_data, scaling_factor, path_to_output):
    locations = gpd.read_file(path_to_locations).to_crs(WGS_84).set_index("id")
    plants = xr.open_dataset(path_to_hydro_data)
    plants = plants.sel(id=(plants.type == "HROR"))
    plant_centroids = gpd.GeoDataFrame(
        crs=WGS_84,
        geometry=list(map(Point, zip(plants.lon, plants.lat))),
        index=plants.id
    )
    location_of_plant = xr.DataArray(
        data=gpd.sjoin(plant_centroids, locations, how="left", op='intersects')["index_right"],
        name="location_id",
        dims=["id"],
        coords={"id": plants["id"]}
    )
    (plants.inflow_MWh.groupby(location_of_plant)
                      .sum(dim="id")
                      .to_dataframe()
                      .reset_index()
                      .pivot(index="time", columns="location_id", values="inflow_MWh")
                      .reindex(columns=locations.index, fill_value=0)
                      .mul(scaling_factor)
                      .rename(columns=lambda col: col.replace(".", "-"))
                      .to_csv(path_to_output, header=True, index=True))


if __name__ == "__main__":
    energy_inflow_per_location(
        path_to_locations=snakemake.input.locations,
        path_to_hydro_data=snakemake.input.hydro,
        scaling_factor=snakemake.params.scaling_factor,
        path_to_output=snakemake.output[0]
    )
