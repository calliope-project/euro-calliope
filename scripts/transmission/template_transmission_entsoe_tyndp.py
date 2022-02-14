import pandas as pd

from eurocalliopelib.template import parametrise_template
from eurocalliopelib import utils


def construct_links(
    path_to_locations, path_to_template, path_to_entsoe_tyndp, scenario, grid, year, ntc_limit,
    scaling_factors, energy_cap_limit, path_to_output
):
    locations = pd.read_csv(path_to_locations, index_col="id")
    tyndp_scenarios = pd.read_excel(path_to_entsoe_tyndp, sheet_name="Line", index_col=0)
    ntcs = _entsoe_ntcs(locations, tyndp_scenarios, scenario, grid, year, ntc_limit)
    link_comment = f"""
# {ntc_limit.title()} net transfer capacity between countries according to the
# following ENTSO-E ten-year network development plan 2020 scenario:
# Scenario: {scenario.title()}
# Case: {grid.title()} Grid
# Year: {year}
# Climate Year: 2007
    """
    parametrise_template(
        path_to_template,
        path_to_output,
        links=ntcs,
        scaling_factors=scaling_factors,
        link_comment=link_comment,
        energy_cap_limit=energy_cap_limit
    )


def _entsoe_ntcs(locations, tyndp_scenarios, scenario, grid, year, ntc_limit):
    """Get Net Transfer Capacities (NTCs) according to the ENTSO-E ten-year network nevelopment plan 2020 scenario dataset"""

    sliced_scenarios = tyndp_scenarios[
        (tyndp_scenarios.Scenario.str.lower() == scenario.lower()) &
        (tyndp_scenarios.Year == year) &
        (tyndp_scenarios["Climate Year"] == 2007)  # not entirely sure what this is, but there is the choice between 1982, 1984, and 2007
    ]
    if grid.lower() == "reference":
        sliced_scenarios = sliced_scenarios[
            sliced_scenarios.Case.str.lower() == "reference grid"
        ]
    elif grid.lower == "expanded":
        sliced_scenarios = sliced_scenarios[
            sliced_scenarios.Case.str.lower().isin(["reference grid", "expanded grid"])
        ]
    # Capacity is in MW
    tyndp_import = _split_links_in_index(sliced_scenarios[sliced_scenarios.Parameter == "Import Capacity"])
    tyndp_export = _split_links_in_index(sliced_scenarios[sliced_scenarios.Parameter == "Export Capacity"])

    # Some NTCs are different depending on whether it is import or export between countries.
    # Here, we take either the minimum or maximum NTC of a link, depending on what a user defines for `ntc_limit`
    average_ntc = pd.concat([tyndp_import, tyndp_export]).abs().agg(ntc_limit, level=["loc_from", "loc_to"])

    # Only keep countries found in the model
    return average_ntc.loc[pd.IndexSlice[locations.index, locations.index]]


def _split_links_in_index(df):
    # Index is LOC[1]-LOC[2], where LOC is is power market
    # (usually a country, but sometimes a subnational group, e.g. 6 in Italy).
    # Here we aggregate to national-level
    df.index = df.index.str.split("-", expand=True)
    df = (
        df
        .rename(index=lambda x: utils.eu_country_code_to_iso3(x[:2]))
        .loc[:, "Value"]
        .groupby(level=[0, 1]).sum()
        .rename_axis(index=["loc_from", "loc_to"])
    )

    # Remove internal transmission for countries with several markets
    return df[df.index.get_level_values(0) != df.index.get_level_values(1)]


if __name__ == "__main__":
    if snakemake.wildcards.resolution == "national":
        construct_links(
            path_to_locations=snakemake.input.locations,
            path_to_template=snakemake.input.template,
            path_to_entsoe_tyndp=snakemake.input.entsoe_tyndp,
            scenario=snakemake.params.scenario,
            grid=snakemake.params.grid,
            year=snakemake.params.year,
            ntc_limit=snakemake.params.ntc_limit,
            scaling_factors=snakemake.params.scaling_factors,
            energy_cap_limit=snakemake.params.energy_cap_limit,
            path_to_output=snakemake.output[0]
        )
