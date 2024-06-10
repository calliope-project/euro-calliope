import numpy as np
import pandas as pd
import pycountry


def scale_to_regional_resolution(df, region_country_mapping, populations, key):
    """
    Create regional vehicle parameters.
    ASSUME all road transport potential is subnationally distributed in proportion to population.
    ASSUME efficiencies of vehicles are equal across regions.
    """

    # Scale charging potential data to population
    if "charging" in key:
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

    # Efficiency data is the same across regions
    if "efficiency" in key:
        regional_df = pd.DataFrame(
            index=df.index,
            data={
                id: df[country_code]
                for id, country_code in region_country_mapping.items()
            },
        ).rename(columns=lambda col_name: col_name.replace(".", "-"))
    return regional_df


def scale_to_national_resolution(df, key):
    df.columns.name = None
    return df


def scale_to_continental_resolution(df, key):
    # Sum all countries' potentials for charging
    if "charging" in key:
        return df.sum(axis=1).to_frame("EUR")
    # Mean of all countries for efficiencies
    if "efficiency" in key:
        return df.mean(axis=1).to_frame("EUR")


def group_vehicle_types(df, vehicle_type_aggregation):
    breakpoint()
    df.index = df.index.set_levels(
        [df.index.levels[0], df.index.levels[1].map(vehicle_type_aggregation)], level=1
    )
    return df.groupby(level=[0, 1]).sum()


def get_efficiency_data_carrier(
    path_to_road_distance: str,
    carrier: str,
    power_scaling_factor: float,
    transport_scaling_factor: float,
    first_year: int,
    final_year: int,
    conversion_factors: dict[str, float],
    country_codes: list[str],
    vehicle_type_aggregation: dict[str, str],
) -> pd.DataFrame:
    """
    Get vehicle efficiency data for the vehicle type corresponding to the carrier, per year
    ASSUME that post 2015 data is the same as 2015 data
    ASSUME that vehicles are grouped according to distance per vehicle type
    """

    def _get_efficiency(x, fuel):
        vehicle_type = x.name[1]
        efficiency = conversion_factors.get(vehicle_type, {}).get(fuel, np.nan)
        return (
            pd.Series(data=[1 / efficiency] * len(x), index=x.index)
            if not pd.isna(efficiency)
            else pd.Series(data=[np.nan] * len(x), index=x.index)
        )

    def _get_distance(vehicle):
        return (
            pd.read_csv(path_to_road_distance)
            .set_index(["vehicle_type", "vehicle_subtype", "year"])
            .xs(slice(first_year, final_year), level="year", drop_level=False)
            .groupby(["vehicle_type", "vehicle_subtype", "country_code", "year"])
            .sum()
            .squeeze()
            .unstack("country_code")
            .loc[:, country_codes]
            .swaplevel(0, 1)
            .loc[[vehicle], :]
        )

    def _aggregate(df, vehicle_type_aggregation):
        # First get the weights of each vehicle category per year
        weights = (
            pd.read_csv(path_to_road_distance)
            .set_index(["vehicle_type", "vehicle_subtype", "year"])
            .xs(slice(first_year, final_year), level="year", drop_level=False)
            .groupby(["vehicle_type", "country_code", "year"])
            .sum()
            .squeeze()
            .unstack("country_code")
            .loc[df.index.get_level_values("vehicle_type").unique(), country_codes]
            .pipe(lambda df: df.div(df.groupby(level="year").sum(), level="year"))
            .pipe(  # Add vehicle_type column
                lambda df: df.assign(
                    vehicle_type_agg=df.index.get_level_values("vehicle_type").map(
                        vehicle_type_aggregation
                    )
                )
            )
            .reset_index()
            .set_index(["year", "vehicle_type", "vehicle_type_agg"])
            .pipe(
                lambda df: df.div(
                    df.groupby(level="vehicle_type_agg").transform("sum"),
                    level="vehicle_type_agg",
                )
            )
        )

        # Then apply the weights to the efficiency data and return
        merged_df = pd.merge(
            df, weights, left_index=True, right_index=True, suffixes=("_df", "_weights")
        )
        for col in df.columns:
            merged_df[col] = merged_df[col + "_df"] * merged_df[col + "_weights"]

        return (
            merged_df.groupby(level=["year", "vehicle_type_agg"])
            .apply(lambda x: x.sum() / x["IRL_weights"].sum())
            .loc[:, country_codes]
        )

    # Return efficiency data for vehicles using the carrier
    carrier_to_vehicle_dict = {
        "electricity": "Battery electric vehicles",
        "diesel": "Diesel oil engine",
    }

    efficiency = (
        _get_distance(carrier_to_vehicle_dict[carrier])
        .apply(_get_efficiency, axis=1, args=(carrier,))
        .droplevel(0)
        .swaplevel(0, 1)
        .pipe(  # Add vehicle_type column
            lambda df: df.assign(
                vehicle_type_agg=df.index.get_level_values("vehicle_type").map(
                    vehicle_type_aggregation
                )
            )
        )
        .reset_index()
        .set_index(["year", "vehicle_type", "vehicle_type_agg"])
    )

    # Aggregate similar vehicle types
    efficiency = _aggregate(efficiency, vehicle_type_aggregation)

    # Missing years are assumed to be the same as the last available year
    if final_year > 2015:
        for year in range(2016, final_year + 1):
            for vehicle_type in efficiency.index.get_level_values(
                "vehicle_type_agg"
            ).unique():
                efficiency.loc[(year, vehicle_type), :] = efficiency.loc[
                    (2015, vehicle_type), :
                ]
    # Return data with correct units
    return efficiency.mul(transport_scaling_factor / power_scaling_factor)


def extract_national_ev_charging_potentials(
    path_to_ev_numbers: str,
    transport_scaling_factor: float,
    first_year: int,
    final_year: int,
    conversion_factors: dict[str, float],
    battery_sizes: dict[str, float],
    country_codes: list[str],
) -> pd.DataFrame:
    # Extract number of EVs per vehicle type
    df_ev_numbers = (
        pd.read_csv(path_to_ev_numbers, index_col=[0, 1, 2, 3, 4])
        .squeeze()
        .droplevel(["vehicle_subtype", "section"])
    )
    assert not df_ev_numbers.isnull().values.any()

    # Compute max. distance travelled per full battery for one EV [in Mio km / vehicle]
    battery_size = pd.DataFrame.from_dict(
        {
            vehicle: battery_sizes[vehicle] / conversion_factors[vehicle]
            for vehicle in battery_sizes
        },
        orient="index",
        columns=["value"],
    )

    # Compute available chargeable distance per vehicle type [in transport scaling unit km]
    df_ev_chargeable_distance = (
        df_ev_numbers.align(battery_size, level="vehicle_type")[1]
        .squeeze()
        .mul(df_ev_numbers)
        .groupby(level=["country_code", "year"])
        .sum()
        .loc[country_codes]
        .unstack("year")
        .mul(transport_scaling_factor)
    )

    if final_year > 2015:
        # ASSUME 2015 data is used for all years after 2015
        df_ev_chargeable_distance = df_ev_chargeable_distance.assign(**{
            str(year): df_ev_chargeable_distance[2015]
            for year in range(2016, final_year + 1)
        })

    df_ev_chargeable_distance.columns = df_ev_chargeable_distance.columns.astype(int)

    return df_ev_chargeable_distance[range(first_year, final_year + 1)].T


def reshape_and_add_suffix(df, suffix):
    return df.T.add_suffix(suffix)


if __name__ == "__main__":
    resolution = snakemake.wildcards.resolution

    power_scaling_factor = snakemake.params.power_scaling_factor
    first_year = snakemake.params.first_year
    final_year = snakemake.params.final_year
    conversion_factors = snakemake.params.conversion_factors
    path_to_output = snakemake.output[0]
    country_codes = [
        pycountry.countries.lookup(c).alpha_3 for c in snakemake.params.countries
    ]
    region_country_mapping = (
        pd.read_csv(snakemake.input.locations, index_col=0)
        .loc[:, "country_code"]
        .to_dict()
    )
    populations = pd.read_csv(snakemake.input.populations, index_col=0)
    battery_sizes = snakemake.params.battery_sizes
    transport_scaling_factor = snakemake.params.transport_scaling_factor
    vehicle_type_aggregation = snakemake.params.vehicle_type_aggregation

    path_to_ev_numbers = snakemake.input.ev_vehicle_number
    path_to_road_distance = snakemake.input.jrc_road_distance

    # Extract efficiency data for electric vehicles
    df_efficiency_ev = get_efficiency_data_carrier(
        path_to_road_distance,
        "electricity",
        power_scaling_factor,
        transport_scaling_factor,
        first_year,
        final_year,
        conversion_factors,
        country_codes,
        vehicle_type_aggregation,
    )

    df_efficiency_ice = get_efficiency_data_carrier(
        path_to_road_distance,
        "diesel",
        power_scaling_factor,
        transport_scaling_factor,
        first_year,
        final_year,
        conversion_factors,
        country_codes,
        vehicle_type_aggregation,
    )

    # Extract national EV charging potentials
    df_charging_potentials = extract_national_ev_charging_potentials(
        path_to_ev_numbers,
        transport_scaling_factor,
        first_year,
        final_year,
        {a: cf["electricity"] for a, cf in conversion_factors.items()},
        battery_sizes,
        country_codes,
    )

    # Add prefix for yaml template
    parameters_evs = {
        "_charging_potential": df_charging_potentials,
        "_efficiency_ev_light": df_efficiency_ev.xs(
            "light_transport", level="vehicle_type_agg"
        ),
        "_efficiency_ev_heavy": df_efficiency_ev.xs(
            "heavy_transport", level="vehicle_type_agg"
        ),
        "_efficiency_ice_light": df_efficiency_ice.xs(
            "light_transport", level="vehicle_type_agg"
        ),
        "_efficiency_ice_heavy": df_efficiency_ice.xs(
            "heavy_transport", level="vehicle_type_agg"
        ),
    }

    # Rescale to desired resolution and add suffix
    # FIXME: rescaling is different for efficiencies and is currently incorrect
    dfs = []
    for key, df in parameters_evs.items():
        if resolution == "continental":
            df = scale_to_continental_resolution(df, key)
        elif resolution == "national":
            df = scale_to_national_resolution(df, key)
        elif resolution in ["regional", "ehighways"]:
            df = scale_to_regional_resolution(
                df,
                region_country_mapping=region_country_mapping,
                populations=populations,
                key=key,
            )
        else:
            raise ValueError("Input resolution is not recognised")
        dfs.append(reshape_and_add_suffix(df, key))

    # Export to csv
    pd.concat(dfs, axis=1).to_csv(path_to_output, index_label=["id"])
