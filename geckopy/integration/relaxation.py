# Copyright 2021 Ginkgo Bioworks

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#     http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Integration layer with pytfa. Pytfa is not installed by default."""

import logging
from typing import Set, TypeVar

import pytfa
import sympy
from pytfa.optim.config import dg_relax_config
from pytfa.optim.constraints import NegativeDeltaG
from pytfa.optim.utils import get_solution_value_for_variables
from pytfa.optim.variables import DeltaGstd, NegSlackVariable, PosSlackVariable
from pytfa.utils.numerics import BIGM_DG
from tqdm import tqdm

import geckopy
from geckopy.experimental.relaxation import Objective_rule, change_constraint


LOGGER = logging.getLogger(__name__)

# this should be a double (AND) constrain not an OR
M = TypeVar("M", geckopy.Model, pytfa.ThermoModel)


def relax_thermo_proteins(
    model: M,
    prot_candidates: Set[str],
    objective_rule: Objective_rule = Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
) -> Set:
    r"""Get one IIS of protein concentrations and the :math:`\Delta G_r` constraints.

    __Inplace__ operation. This combines both `pytfa.optim.relax_dgo` and
    :py:func:`geckopy.experimental.relaxation.get_upper_relaxation` and avoids
    copying the model since that fails for a `pytfa.ThermoModel` with proteins.

    Parameters
    ----------
    original_model: cobra.Model
    elastic_candidates: list[str]
    objective_rule: Objective_rule
        The IIS is selected by minimizing an objective as defined in
        :class:`Objective_rule`.

    Returns
    -------
    iis: Set
        variables that relax the model to make it feashible.
    status: str
    """
    dg_relax_config(model)
    # 1. add elastic thermodymanic variables
    my_dgo = model.get_variables_of_type(DeltaGstd)
    my_neg_dg = model.get_constraints_of_type(NegativeDeltaG)
    objective_symbols = []
    model.solver.update()
    for this_neg_dg in tqdm(my_neg_dg, desc="adding thermo slacks"):
        if this_neg_dg.id not in my_dgo:
            continue
        # Create the negative and positive slack variables
        # We can't queue them because they will be in an expression to declare
        # the constraint
        neg_slack = model.add_variable(
            NegSlackVariable, this_neg_dg.reaction, lb=0, ub=BIGM_DG, queue=False
        )
        pos_slack = model.add_variable(
            PosSlackVariable, this_neg_dg.reaction, lb=0, ub=BIGM_DG, queue=False
        )

        # Create the new constraint by adding the slack variables to the
        # negative delta G constraint (from the initial cobra_model)
        new_expr = this_neg_dg.constraint.expression
        new_expr += pos_slack - neg_slack

        this_neg_dg.change_expr(new_expr)

        # Update the objective with the new variables
        objective_symbols += [neg_slack.variable, pos_slack.variable]
    # 2. add elastic upper protein variables
    for e_id in tqdm(prot_candidates, desc="adding protein slacks"):
        e = model.constraints[e_id]
        if objective_rule == Objective_rule.MIN_MILP_COUNT:
            v = model.problem.Variable(f"V_{e_id}", type="binary")
            vl = model.problem.Variable(f"VL_{e_id}", lb=0, ub=1000)
            new_expr = e.expression + v * 1000 - vl
        else:
            v = model.problem.Variable(f"V_{e_id}", lb=0, ub=1000)
            new_expr = e.expression + v
        change_constraint(e, new_expr, model.solver)
        objective_symbols.append(v)

    # Change the objective to minimize slack
    if objective_rule == Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE:
        model.objective = model.problem.Objective(
            sympy.Add(*(objective_symbols + [-model.objective.expression])),
            direction="min",
        )
    else:
        model.objective = model.problem.Objective(
            sympy.Add(*(objective_symbols)),
            direction="min",
        )
    model.repair()
    relaxation = model.optimize()
    status = model.solver.status
    my_neg_slacks = model.get_variables_of_type(NegSlackVariable)
    my_pos_slacks = model.get_variables_of_type(PosSlackVariable)
    # filter by all three elastic variable prefixes used
    neg_series = get_solution_value_for_variables(relaxation, my_neg_slacks)
    pos_series = get_solution_value_for_variables(relaxation, my_pos_slacks)
    return (
        {
            v.name[2:]
            for v in objective_symbols
            if abs(v.primal) > 0 and (v.name.startswith("V_"))
        }
        | set(neg_series[neg_series.abs() > 0].index)
        | set(pos_series[pos_series.abs() > 0].index)
        if status == "optimal"
        else set(),
        status,
    )
