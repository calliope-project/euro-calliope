import pandas as pd
from eurocalliopelib.template import parametrise_template

if __name__ == "__main__":
    biofuel_cost = float(pd.read_csv(snakemake.input.biofuel_cost).columns[0])
    locations = pd.read_csv(snakemake.input.locations, index_col=0)

    parametrise_template(
        snakemake.input.template,
        snakemake.output[0],
        biofuel_cost=biofuel_cost,
        scaling_factors=snakemake.params.scaling_factors,
        locations=locations,
    )
