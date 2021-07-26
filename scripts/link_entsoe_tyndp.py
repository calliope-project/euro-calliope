

import pandas as pd
import jinja2

from eurocalliopelib import filters, utils

TEMPLATE = """
# {{ entsoe_ntc_limit.title() }} net transfer capacity between countries according to the
# following ENTSO-E ten-year network development plan 2020 scenario:
# Scenario: {{ entsoe_scenario.title() }}
# Case: {{ entsoe_grid.title() }} Grid
# Year: {{ entsoe_year }}
# Climate Year: 2007

links:
    {% for link in ntcs.index %}
    {{ link[0] }},{{ link[1] }}:
        techs.ac_transmission.constraints:
            energy_cap_{{ energy_cap_limit }}: {{ ntcs.loc[link] * scaling_factor }}  # {{ (1 / scaling_factor) | unit("MW") }}
    {% endfor %}
"""

def link_tyndp(
    path_to_locations, path_to_entsoe_tyndp, entsoe_scenario, entsoe_grid, entsoe_year, entsoe_ntc_limit,
    scaling_factor, energy_cap_limit, path_to_result
):
    locations = pd.read_csv(path_to_locations, index_col="id")
    ntcs = _entsoe_ntcs(locations, path_to_entsoe_tyndp, entsoe_scenario, entsoe_grid, entsoe_year, entsoe_ntc_limit)

    env = jinja2.Environment()
    env.filters['unit'] = filters.unit
    links = env.from_string(TEMPLATE).render(
        ntcs=ntcs,
        scaling_factor=scaling_factor,
        energy_cap_limit=energy_cap_limit,
        entsoe_scenario=entsoe_scenario,
        entsoe_grid=entsoe_grid,
        entsoe_year=entsoe_year,
        entsoe_ntc_limit=entsoe_ntc_limit,
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(links)


def _entsoe_ntcs(locations, path_to_entsoe_tyndp, entsoe_scenario, entsoe_grid, entsoe_year, entsoe_ntc_limit):
    """Get Net Transfer Capacities (NTCs) according to the ENTSO-E ten-year network nevelopment plan 2020 scenario dataset"""
    tyndp_scenarios = pd.read_excel(path_to_entsoe_tyndp, sheet_name="Line", index_col=0)
    tyndp_scenarios = tyndp_scenarios[
        (tyndp_scenarios.Scenario.str.lower() == entsoe_scenario.lower()) &
        (tyndp_scenarios.Case.str.lower() == f"{entsoe_grid} Grid".lower()) &
        (tyndp_scenarios.Year == entsoe_year) &
        (tyndp_scenarios["Climate Year"] == 2007)  # not entirely sure what this is, but there is the choice between 1982, 1984, and 2007
    ]
    # Capacity is in MW
    tyndp_import = _split_links_in_index(tyndp_scenarios[tyndp_scenarios.Parameter == "Import Capacity"])
    tyndp_export = _split_links_in_index(tyndp_scenarios[tyndp_scenarios.Parameter == "Export Capacity"])

    # Some NTCs are different depending on whether it is import or export between countries.
    # Here, we take either the minimum or maximum NTC of a link, depending on what a user defines for `entsoe_ntc_limit`
    average_ntc = getattr(pd.concat([tyndp_import, tyndp_export]).abs(), entsoe_ntc_limit)(level=["loc_from", "loc_to"])

    # Only keep countries found in the model
    return average_ntc[
        average_ntc.index.get_level_values(0).isin(locations.index) &
        average_ntc.index.get_level_values(1).isin(locations.index)
    ]


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
        link_tyndp(
            path_to_locations=snakemake.input.units,
            path_to_entsoe_tyndp=snakemake.input.entsoe_tyndp,
            entsoe_scenario=snakemake.params.entsoe_scenario,
            entsoe_grid=snakemake.params.entsoe_grid,
            entsoe_year=snakemake.params.entsoe_year,
            entsoe_ntc_limit=snakemake.params.entsoe_ntc_limit,
            scaling_factor=snakemake.params.scaling_factor,
            energy_cap_limit=snakemake.params.energy_cap_limit,
            path_to_result=snakemake.output[0]
        )
