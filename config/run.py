import sys

import calliope
import numpy as np
import pyomo.core as po
from calliope.backend.pyomo.util import (
    get_param,
    invalid,
    split_comma_list,
)

"""
This file allows to run a Calliope model with custom constraints

To run it, use the following command:
    python run.py path_to_model scenario_string path_to_output

    where:
    - path_to_model is the path to the Calliope model to run,
    - scenario_string is the scenario to run,
    - and path_to_output is the path to the output .nc file.

Make sure to have a Calliope environment which can run the model, e.g. the test environment.
"""

args = sys.argv[1:]
path_to_model = args[0]
scenario_string = args[1]
path_to_output = args[2]


def build_model(path_to_model, scenario, path_to_temp_output):
    # Build model without custom constraints. Does not run the model
    calliope.set_log_verbosity(
        "info", include_solver_output=True, capture_warnings=True
    )
    model = calliope.Model(path_to_model, scenario=scenario)

    model._model_data.attrs["scenario"] = scenario

    model.to_netcdf(path_to_temp_output)


def run_model(path_to_temp_model, path_to_output):
    # Run model with custom constraints
    calliope.set_log_verbosity(
        "info", include_solver_output=True, capture_warnings=True
    )
    model = calliope.read_netcdf(path_to_temp_model)
    model.run(build_only=True)
    add_eurocalliope_constraints(model)
    new_model = model.backend.rerun()
    new_model.to_netcdf(path_to_output)


def add_eurocalliope_constraints(model):
    backend_model = model._backend_model
    if "energy_cap_max_time_varying" in model._model_data.data_vars:
        print("Building production_max_time_varying constraint")
        add_production_max_time_varying_constraint(model, backend_model)
    if any(
        var.startswith("demand_shape_per_month") for var in model._model_data.data_vars
    ):
        print("Building demand_shape_per_month constraint")
        add_carrier_prod_per_month_constraints(model, backend_model)


def equalizer(lhs, rhs, sign):
    if sign == "max":
        return lhs <= rhs
    elif sign == "min":
        return lhs >= rhs
    elif sign == "equals":
        return lhs == rhs
    else:
        raise ValueError(f"Invalid sign: {sign}")


def add_production_max_time_varying_constraint(model, backend_model):
    def _carrier_production_max_time_varying_constraint_rule(
        backend_model, loc_tech, timestep
    ):
        """
        Set maximum carrier production for technologies with time varying maximum capacity
        """
        energy_cap_max = backend_model.energy_cap_max_time_varying[loc_tech, timestep]
        if invalid(energy_cap_max):
            return po.Constraint.Skip
        model_data_dict = backend_model.__calliope_model_data["data"]
        timestep_resolution = backend_model.timestep_resolution[timestep]
        loc_tech_carriers_out = split_comma_list(
            model_data_dict["lookup_loc_techs_conversion"]["out", loc_tech]
        )

        carrier_prod = sum(
            backend_model.carrier_prod[loc_tech_carrier, timestep]
            for loc_tech_carrier in loc_tech_carriers_out
        )
        return carrier_prod <= (
            backend_model.energy_cap[loc_tech] * timestep_resolution * energy_cap_max
        )

    backend_model.loc_tech_carrier_production_max_time_varying_constraint = po.Set(
        initialize=[
            loc_tech
            for loc_tech in backend_model.loc_techs_conversion
            if model.inputs.energy_cap_max_time_varying.loc[{"loc_techs": loc_tech}]
            .notnull()
            .all()
        ],
        ordered=True,
    )
    model.backend.add_constraint(
        "carrier_production_max_time_varying_constraint",
        ["loc_tech_carrier_production_max_time_varying_constraint", "timesteps"],
        _carrier_production_max_time_varying_constraint_rule,
    )


def add_carrier_prod_per_month_constraints(model, backend_model):
    def _carrier_prod_per_month_constraint_rule_generator(sense):
        def __carrier_prod_per_month_constraint_rule(backend_model, loc_tech, month):
            """
            Set the min/max amount of carrier consumption (relative to annual consumption)
            for a specific loc tech that must take place in a given calender month in the model
            """
            model_data_dict = backend_model.__calliope_model_data
            loc_tech_carrier = model_data_dict["data"]["lookup_loc_techs_conversion"][
                ("out", loc_tech)
            ]

            prod = backend_model.carrier_prod
            prod_total = sum(
                prod[loc_tech_carrier, timestep] for timestep in backend_model.timesteps
            )
            prod_month = sum(
                prod[loc_tech_carrier, timestep]
                for timestep in backend_model.timesteps
                if backend_model.month_numbers[timestep].value == month
            )
            if "timesteps" in [
                i.name
                for i in getattr(
                    backend_model, f"carrier_prod_per_month_{sense}_time_varying"
                )._index.subsets()
            ]:
                prod_fraction = sum(
                    get_param(
                        backend_model,
                        f"carrier_prod_per_month_{sense}_time_varying",
                        (loc_tech, timestep),
                    )
                    * backend_model.timestep_resolution[timestep]
                    for timestep in backend_model.timesteps
                    if backend_model.month_numbers[timestep].value == month
                )
            else:
                prod_fraction = get_param(
                    backend_model, f"carrier_prod_per_month_{sense}", (loc_tech)
                )

            return equalizer(prod_month, prod_total * prod_fraction, sense)

        return __carrier_prod_per_month_constraint_rule

    backend_model.months = po.Set(
        initialize=np.unique(model._model_data.timesteps.dt.month.values), ordered=True
    )
    month_numbers = model._model_data.timesteps.dt.month.to_series()
    month_numbers.index = month_numbers.index.strftime("%Y-%m-%d %H:%M")

    backend_model.month_numbers = po.Param(
        backend_model.timesteps,
        initialize=month_numbers.to_dict(),
        mutable=True,
        within=po.Reals,
    )
    backend_model.__calliope_datetime_data.add(("data_vars", "month_numbers"))

    for sense in ["min", "max", "equals"]:
        if (
            f"carrier_prod_per_month_{sense}_time_varying"
            in model._model_data.data_vars
        ):
            setattr(
                backend_model,
                f"loc_techs_carrier_prod_per_month_{sense}",
                po.Set(
                    initialize=[
                        loc_tech
                        for loc_tech in backend_model.loc_techs
                        if (
                            model._model_data[
                                f"carrier_prod_per_month_{sense}_time_varying"
                            ]
                            .loc[{"loc_techs": loc_tech}]
                            .notnull()
                            .all()
                        )
                    ],
                    ordered=True,
                ),
            )
            model.backend.add_constraint(
                f"carrier_prod_per_month_{sense}_constraint",
                [f"loc_techs_carrier_prod_per_month_{sense}", "months"],
                _carrier_prod_per_month_constraint_rule_generator(sense),
            )


build_model(path_to_model, scenario_string, path_to_output)
run_model(path_to_output, path_to_output)

model = calliope.read_netcdf(path_to_output)
model.plot.timeseries(subset={"locs": ["IRL"]})
