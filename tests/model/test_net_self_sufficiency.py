import calliope
import pytest
import pandas as pd

DEMAND_TECH = "demand_elec"
EPSILON_IMPORTS = 0.001 # 0.1 %
NET_SELF_SUFFICIENCY_LEVELS = [100, 85]


@pytest.fixture
def self_sufficient_model(path_to_test_model, net_self_sufficiency_level):
    return calliope.Model(
        path_to_test_model,
        scenario=f"regional-net-self-sufficiency-{net_self_sufficiency_level}-percent"
    )


@pytest.fixture(params=NET_SELF_SUFFICIENCY_LEVELS)
def net_self_sufficiency_level(request):
    return request.param


@pytest.fixture
def transmission_loc_techs_per_self_sufficient_group(self_sufficient_model):
    return [
        getattr(self_sufficient_model._model_data, f"group_constraint_loc_techs_{group}")
        for group in self_sufficient_model._model_data.group_names_net_import_share_max.values
    ]


@pytest.fixture
def net_imports_per_self_sufficient_group(self_sufficient_model, transmission_loc_techs_per_self_sufficient_group, scaling_factors):
    imports = [
        (self_sufficient_model
         .results
         .carrier_prod
         .sel(loc_tech_carriers_prod=[f"{loc_tech.item()}::electricity" for loc_tech in transmission_loc_techs])
         .sum()
         .item()) / scaling_factors["power"]
        for transmission_loc_techs in transmission_loc_techs_per_self_sufficient_group
    ]
    exports = [
        (self_sufficient_model
         .results
         .carrier_con
         .sel(loc_tech_carriers_con=[f"{loc_tech.item()}::electricity" for loc_tech in transmission_loc_techs])
         .sum()
         .item()) / scaling_factors["power"]
        for transmission_loc_techs in transmission_loc_techs_per_self_sufficient_group
    ]
    return [i + e for i, e in zip(imports, exports)]


@pytest.fixture
def demand_per_self_sufficient_group(self_sufficient_model, transmission_loc_techs_per_self_sufficient_group):
    carrier_con = self_sufficient_model.get_formatted_array("carrier_con")
    locs_per_autarkic_group = [
        list(set([loc_tech.item().split("::")[0] for loc_tech in transmission_loc_techs]))
        for transmission_loc_techs in transmission_loc_techs_per_self_sufficient_group
    ]
    return [
        carrier_con.sel(techs=DEMAND_TECH, locs=locs).sum().item()
        for locs in locs_per_autarkic_group
    ]


def test_net_imports(net_imports_per_self_sufficient_group, demand_per_self_sufficient_group, rel_net_import_threshold):
    rel_net_imports_per_self_sufficient_group = pd.Series([
        i / -d
        for i, d in zip(net_imports_per_self_sufficient_group, demand_per_self_sufficient_group)
    ])
    assert (rel_net_imports_per_self_sufficient_group <= (rel_net_import_threshold + EPSILON_IMPORTS)).all()
