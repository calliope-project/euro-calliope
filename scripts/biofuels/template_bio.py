import pandas as pd
from eurocalliopelib.template import parametrise_template


def construct_techs_and_locations(
    path_to_template,
    path_to_output,
    path_to_locations,
    path_to_biofuel_costs,
    biofuel_efficiency,
    scaling_factors,
):
    with open(path_to_biofuel_costs) as f_biofuel_costs:
        biofuel_fuel_cost = float(f_biofuel_costs.readline())
    locations = pd.read_csv(path_to_locations, index_col=0)

    return parametrise_template(
        path_to_template,
        path_to_output,
        biofuel_fuel_cost=biofuel_fuel_cost,
        biofuel_efficiency=biofuel_efficiency,
        scaling_factors=scaling_factors,
        locations=locations,
    )


if __name__ == "__main__":
    construct_techs_and_locations(
        path_to_template=snakemake.input.template,
        path_to_locations=snakemake.input.locations,
        path_to_biofuel_costs=snakemake.input.biofuel_cost,
        biofuel_efficiency=snakemake.params.biofuel_efficiency,
        scaling_factors=snakemake.params.scaling_factors,
        path_to_output=snakemake.output[0],
    )
