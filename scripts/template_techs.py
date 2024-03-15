"""
Creates Calliope tech files, including allocation to locations and overrides, when there are no special inputs being supplied.
"""

import pandas as pd
from eurocalliopelib.template import parametrise_template


def construct_techs_and_locations(
    path_to_template, path_to_output, path_to_locations, params
):
    locations = pd.read_csv(path_to_locations, index_col=0)
    parametrise_template(
        path_to_template, path_to_output, locations=locations, **params
    )


if __name__ == "__main__":
    construct_techs_and_locations(
        path_to_template=snakemake.input.template,
        path_to_locations=snakemake.input.locations,
        params=snakemake.params,
        path_to_output=snakemake.output[0],
    )
