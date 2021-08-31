import jinja2
import pandas as pd

TEMPLATE = """
overrides:
    {%- for level in levels %}
    {%- for scale, units in groups.items() %}
    {{ scale }}-{{ level }}-percent:
        {%- for group_name, sublocations in units.items() %}
        group_constraints.self_sufficiency_{{ level }}_percent_{{ group_name }}:
            locs:
            {%- for sublocation in sublocations %}
            - {{ sublocation }}
            {%- endfor %}
            net_import_share_max:
                electricity: {{ (1 - level / 100) | round(2) }}
        {%- endfor %}
    {%- endfor %}
    {%- endfor %}
"""
COUNTRY_CODE_COLUMN = "country_code"


def net_self_sufficiency(path_to_units, levels_in_percent, connected_regions, path_to_result):
    """Create overrides for net self-sufficicency on different geographic scales and with different levels."""
    units = pd.read_csv(path_to_units, index_col=0).rename(index=lambda idx: idx.replace(".", "-"))

    groups = {
        "regional-net-self-sufficiency": _regional_self_sufficiency(units),
        "national-net-self-sufficiency": _national_self_sufficiency(units),
        "continental-net-self-sufficiency": _continental_self_sufficiency(units)
    }
    groups["regional-net-self-sufficiency"] = _connect_regions(
        groups["regional-net-self-sufficiency"],
        connected_regions
    )

    restrictions = jinja2.Template(TEMPLATE).render(
        groups=groups,
        levels=levels_in_percent
    )
    with open(path_to_result, "w") as result_file:
        result_file.write(restrictions)


def _regional_self_sufficiency(units):
    # assume units are regions or larger
    return {unit_id: [unit_id] for unit_id, unit in units.iterrows()}


def _national_self_sufficiency(units):
    return {
        country_code: list(countries.index.values)
        for country_code, countries in units.groupby(COUNTRY_CODE_COLUMN)
    }


def _continental_self_sufficiency(units):
    return {
        "EUR": list(units.index.values)
    }


def _connect_regions(groups, connected_regions):
    if all([region1 in groups.keys() and region2 in groups.keys()
            for region1, region2 in connected_regions]):
        # config is valid, and resolution is regional. apply config
        for region1, region2 in connected_regions:
            groups[region1] = [region1, region2]
            del groups[region2]
        return groups
    elif all([region1 not in groups.keys() and region2 not in groups.keys()
              for region1, region2 in connected_regions]):
        # config is valid, but resolution is not regional. do nothing
        return groups
    else:
        raise ValueError("Config of connected regions is invalid.")


if __name__ == "__main__":
    net_self_sufficiency(
        path_to_units=snakemake.input.units,
        levels_in_percent=snakemake.params.levels,
        connected_regions=snakemake.params.connected_regions,
        path_to_result=snakemake.output[0]
    )
