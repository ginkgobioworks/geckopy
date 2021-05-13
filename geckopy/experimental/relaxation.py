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

"""Relax experimental constraints.

Usually, the experimental constraints produce an infeashible model. The
relaxation methods aim to remove the smallest subset of experimental measurements
to allow growth.
"""

import logging
import warnings
from enum import Enum
from math import isnan
from typing import Dict, List, Optional, Set, Tuple

import cobra
import optlang
import pandas as pd
import sympy

from geckopy.model import Model


__all__ = [
    "apply_proteomics_relaxation",
    "apply_proteomics_elastic_relaxation",
    "elastic_upper_relaxation",
    "relax_proteomics_greedy",
    "get_upper_relaxation",
    "Objective_rule",
]
LOGGER = logging.getLogger(__name__)


class Objective_rule(Enum):
    r"""Objective to minimize for relaxation.

    - Objective_rule.MIN_ELASTIC_SUM:
    :math:`\sum_{v \in \text{elastic vars}} v_{flux}` (LP).
    - Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE:
    :math:`\sum_{e \in \text{elastic vars}} + \text{prev objective}` (LP).
    - Objective_rule.MIN_MILP_COUNT:
    :math:`\sum_i^{N} e_i` where e is a binary variable (MILP).
    """

    MIN_ELASTIC_SUM = 1
    MIN_MILP_COUNT = 2
    MIN_ELASTIC_SUM_OBJECTIVE = 3


def apply_proteomics_relaxation(
    original_model: Model,
    objective_rule: Objective_rule = Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
) -> Tuple[Model, Set]:
    """Relax the problem by relaxing the protein concentration constraints.

    The relaxed problems will contain elastic variables, returning the model a
    non-unique infeashible contraint set.

    Parameters
    ----------
    original_model: geckopy.Model
        Geckopy model. It won't be modified but copied.
    objective_rule: Objective_rule
        The IIS is selected by minimizing an objective as defined in
        :class:`Objective_rule`.

    Returns
    -------
    tuple: (geckopy.Model, set)
        copy of the model with the relaxed variables applied and the sets

    """
    model = original_model.copy()
    iis, _ = get_upper_relaxation(
        model,
        [
            prot.id
            for prot in model.proteins
            if isinstance(prot.concentration, float) and not isnan(prot.concentration)
        ],
        objective_rule,
    )
    # second relaxation to generate the model only with the required new
    # variables and with the modified objective
    return model, iis


def apply_proteomics_elastic_relaxation(
    original_model: Model,
    objective_rule: Objective_rule = Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
) -> Tuple[Model, Set]:
    """Relax the problem by relaxing the protein concentration constraints.

    The relaxed problems will be determined via Elastic filtering, returning
    the model and a irreducibly inconsistent set of functional constraints (
    [Chinnek and Dravnieks, 1990]
    (https://pubsonline.informs.org/doi/abs/10.1287/ijoc.3.2.157)).

    Parameters
    ----------
    original_model: geckopy.Model
        Geckopy model. It won't be modified but copied.
    objective_rule: Objective_rule
        The IIS is selected by minimizing an objective as defined in
        :class:`Objective_rule`.

    Returns
    -------
    tuple: (geckopy.Model, set)
        copy of the model with the relaxed variables applied and the sets
    """
    model = original_model.copy()
    # model is inspescted for IIS
    elastics = elastic_upper_relaxation(
        model,
        [
            prot.id
            for prot in model.proteins
            if isinstance(prot.concentration, float) and not isnan(prot.concentration)
        ],
        objective_rule,
    )
    # model is modified in place given the elastic candidates found
    iis, _ = get_upper_relaxation(model, elastics, objective_rule)
    return model, iis


def change_constraint(
    this_constr: optlang.Constraint,
    new_expr: optlang.Constraint,
    model,
    sloppy: bool = False,
    lb: Optional[float] = None,
    ub: Optional[float] = None,
):
    """Change a constraint expression in the solver."""
    lb = this_constr.lb if lb is None else lb
    ub = this_constr.ub if ub is None else ub
    name = this_constr.name
    # Remove former constraint to override it
    model.remove(this_constr)
    new_cons = model.interface.Constraint(name=name, expression=new_expr, ub=ub, lb=lb)
    # Add the new variant
    model.add(new_cons, sloppy=sloppy)


def get_upper_relaxation(
    model: cobra.Model,
    candidates: List[str],
    objective_rule: Objective_rule = Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
) -> Tuple[Set, str]:
    """Get one IIS of upper bounds of the `candidates`."""
    objective_vars = []
    for e_id in candidates:
        e = model.constraints[e_id]
        if objective_rule == Objective_rule.MIN_MILP_COUNT:
            v = model.problem.Variable(f"V_{e_id}", type="binary")
            vl = model.problem.Variable(f"VL_{e_id}", lb=0, ub=1000)
            new_expr = e.expression + v * 1000 - vl
        else:
            v = model.problem.Variable(f"V_{e_id}", lb=0, ub=1000)
            new_expr = e.expression + v
        change_constraint(e, new_expr, model.solver)
        objective_vars.append(v)
    if objective_rule == Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE:
        model.objective = model.problem.Objective(
            sympy.Add(*(objective_vars + [-model.objective.expression])),
            direction="min",
        )
    else:
        model.objective = model.problem.Objective(
            sympy.Add(*(objective_vars)),
            direction="min",
        )
    model.slim_optimize()
    status = model.solver.status
    # CPLEX may just error out if a var of an infeasible problem is accessed
    return {
        v.name[2:]
        for v in objective_vars
        if abs(v.primal) > 0 and v.name.startswith("V_")
    } if status == "optimal" else set(), status


def elastic_upper_relaxation(
    original_model: cobra.Model,
    elastic_candidates: List[str],
    objective_rule: Objective_rule = Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
) -> Set:
    r"""Convert constrains to elastic constraints until the problem is feashible.

    It assumes that the elastic candidates are all subject to a <= constraint.

    Based on Brown and Graves, 1975:

        - Q <- original problem

        - E <- set of variables, candidates to be relaxed

        - V <- set of variables that relax the problem

        - relax(:math:`e_i`) <- Add constraint :math:`e - v`

        - Z <- :math:`\sum_{v \in V} v`

    IIS <- irreducibly infeashible set of variables

    1. relax(:math:`e_i`) for each :math:`e_i` in :math:`E`
    2. min :math:`Z`, s.t. :math:`E, Q, v \ge 0 \forall V`
    3. Get :math:`R = v \in V \lt 0`
    4. If infeashible:
        STOP
    5. Else:
        5.1. :math:`E = E - R`

        5.2. :math:`IIS = IIS \cup R`

        Go to 1
    6. return IIS

    Parameters
    ----------
    original_model: cobra.Model
    elastic_candidates: list[str]
    objective_rule: Objective_rule
        The IIS is selected by minimizing an objective as defined in
        :class:`Objective_rule`.
    """
    status = "optimal"
    current_candidates = set(elastic_candidates)
    iss = set()
    n_last_iss = -1
    while status == "optimal" and n_last_iss != len(iss):
        n_last_iss = len(iss)
        current_candidates -= iss
        n_iss, status = get_upper_relaxation(
            original_model.copy(), current_candidates, objective_rule
        )
        # CPLEX may just error out if a var of an infeasible problem is accessed
        iss |= n_iss
    return iss


def top_shadow_prices(
    solution: cobra.Solution, top: int = 1, protein_set: Optional[List[str]] = None
) -> pd.DataFrame:
    """Rank them from most to least sensitive in the model reactions.

    Parameters
    ----------
    solution: cobra.Solution
        The usual Solution object returned by model.optimize().
    top: int
        The number of metabolites to be returned.
    protein_set: Optional[List[str]]
        If a list of ids is provided, the search will be applied to only these proteins.

    Returns
    -------
    shadow_pr: pd.Series
        Top shadow prices, ranked.
    """
    protein_set = (
        protein_set if protein_set is not None else solution.shadow_prices.index
    )
    shadow_pr = solution.shadow_prices
    return shadow_pr[protein_set].sort_values()[:top]


def relax_proteomics_greedy(
    model: Model, minimal_growth: float, protein_set: Optional[List[str]] = None
) -> Tuple[Dict, List[Dict]]:
    """Remove proteomics measurements with a set that enables the model to grow.

    Proteins are removed from the set iteratively based on sensitivity analysis
    (shadow prices).

    Adapted from
    https://github.com/DD-DeCaF/simulations/blob/devel/src/simulations/modeling/driven.py

    Parameters
    ----------
    model: cobra.Model
        The enzyme-constrained model.
    minimal_growth_rate: float
        Minimal growth rate to enforce.
    protein_set: Optional[List[str]]
        If a list of ids is provided, the search will be applied to only these proteins.

    Returns
    -------
    growth_rate: dict
        New growth rate (will change if the model couldn't grow at the inputted value).
    proteomics: list(dict)
        Filtered list of proteomics.

    """
    solution, prot_solution = model.optimize()
    new_growth_rate = (
        solution.objective_value
        if solution.objective_value and not isnan(solution.objective_value)
        else 0
    )

    # while the model cannot grow to the desired level, remove the protein with
    # the highest shadow price:
    prots_to_remove = []
    n_prots_with_concentration = (
        len(
            [
                prot
                for prot in model.proteins
                if prot.concentration and not isnan(prot.concentration)
            ]
        )
        if protein_set is None
        else len(protein_set)
    )
    while new_growth_rate < minimal_growth and n_prots_with_concentration > 0:
        # get most influential protein in model:
        top_protein = top_shadow_prices(prot_solution, protein_set=protein_set)
        top_protein = model.proteins.get_by_id(top_protein.index[0])

        # update data: append protein to list, remove from current dataframe and
        # increase the corresponding upper bound to +1000:
        prots_to_remove.append(top_protein)
        LOGGER.debug(f"Removed {top_protein.id} concentration")
        top_protein.add_concentration(None)
        if protein_set:
            protein_set.remove(top_protein.id)
        n_prots_with_concentration -= 1

        # re-compute solution:
        solution, prot_solution = model.optimize()
        new_growth_rate = (
            solution.objective_value
            if solution.objective_value and not isnan(solution.objective_value)
            else 0
        )

    # update growth rate if optimization was not successful:
    if new_growth_rate < minimal_growth:
        message = (
            "Minimal growth was not reached! "
            f"Final growth of the model: {new_growth_rate}"
        )
        LOGGER.warn(message)
        warnings.warn(message)

    return new_growth_rate, prots_to_remove
