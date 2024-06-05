import pandas as pd
import pycountry


def scale_to_regional_resolution(df, region_country_mapping, populations):
    """
    ASSUME data is subnationally distributed in proportion to population.
    """
    df_population_share = (
        populations.loc[:, "population_sum"]
        .reindex(region_country_mapping.keys())
        .groupby(by=region_country_mapping)
        .transform(lambda df: df / df.sum())
    )

    regional_df = (
        pd.DataFrame(
            index=df.index,
            data={
                id: df[country_code]
                for id, country_code in region_country_mapping.items()
            },
        )
        .mul(df_population_share)
        .rename(columns=lambda col_name: col_name.replace(".", "-"))
    )
    pd.testing.assert_series_equal(regional_df.sum(axis=1), df.sum(axis=1))
    return regional_df


def scale_to_national_resolution(df):
    df.columns.name = None
    return df


def scale_to_continental_resolution(df):
    return df.sum(axis=1).to_frame("EUR")


def rescale_to_resolution(df, resolution):
    if resolution in ["regional", "ehighways"]:
        return scale_to_regional_resolution(
            df_annual_road_distance_demand, region_country_mapping, populations
        )
    elif resolution == "national":
        return scale_to_national_resolution(df_annual_road_distance_demand)
    elif resolution == "continental":
        return scale_to_continental_resolution(df_annual_road_distance_demand)


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution
    path_to_annual_road_distance = snakemake.input.annual_road_distance
    first_year = snakemake.params.first_year
    final_year = snakemake.params.final_year
    path_to_output = snakemake.output[0]
    vehicle_aggregation = snakemake.params.vehicle_aggregation

    country_codes = [
        pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
    ]
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )

    populations = pd.read_csv(snakemake.input.populations, index_col=0)

    # Extract annual road distance demand data
    df_annual_road_distance_demand = (
        pd.read_csv(path_to_annual_road_distance, index_col=[0, 1, 2])
        .squeeze()
        .rename(vehicle_aggregation, axis=0, level="vehicle_type")
        .groupby(level=["vehicle_type", "country_code", "year"])
        .sum()
        .xs(slice(first_year, final_year), level="year", drop_level=False)
        .unstack("country_code")
        .loc[:, country_codes]
    )

    # Rescale annual demand data to resolution and reshape format
    df_annual_road_distance_demand = (
        rescale_to_resolution(df_annual_road_distance_demand, resolution)
        .swaplevel()
        .sort_index(level=0)
        .T
    )

    # Flatten for YAML output
    df_annual_road_distance_demand.columns = [
        f"{col[0]}_{col[1]}" for col in df_annual_road_distance_demand.columns
    ]
    df_annual_road_distance_demand.index.name = "id"

    # Multiply by transport unit and export to CSV
    df_annual_road_distance_demand.mul(snakemake.params.transport_factor).to_csv(
        path_to_output
    )
