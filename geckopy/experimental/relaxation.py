"""Relax experimental constraints.

Usually, the experimental constraints produce an infeashible model. The
relaxation methods aim to remove the smallest subset of experimental measurements
to allow growth.
"""

from math import isnan
from typing import Dict, List, Set, Tuple

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
]


def apply_proteomics_relaxation(
    original_model: Model, min_objective=0.0
) -> Tuple[Model, Set]:
    """Relax the problem by relaxing the protein concentration constraints.

    The relaxed problems will contain elastic variables, returning the model a
    non-minimal infeashible contraint set.

    Parameters
    ----------
    original_model: geckopy.Model
        Geckopy model. It won't be modified but copied.
    min_growth: float
        Minimal current objective value to be reached.

    Returns
    -------
    tuple: (geckopy.Model, set)
        copy of the model with the relaxed variables applied and the sets

    """
    model = original_model.copy()
    _, elastics = apply_upper_relaxation(
        model,
        [
            prot.id
            for prot in model.proteins
            if isinstance(prot.concentration, float) and not isnan(prot.concentration)
        ],
    )
    # second relaxation to generate the model only with the required new
    # variables and with the modified objective
    return apply_upper_relaxation(model, elastics)


def apply_proteomics_elastic_relaxation(
    original_model: Model, min_objective=0.0
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
    min_growth: float
        Minimal current objective value to be reached.

    Returns
    -------
    tuple: (geckopy.Model, set)
        copy of the model with the relaxed variables applied and the sets
    """
    model = original_model.copy()
    elastics = elastic_upper_relaxation(
        model,
        [
            prot.id
            for prot in model.proteins
            if isinstance(prot.concentration, float) and not isnan(prot.concentration)
        ],
    )
    return apply_upper_relaxation(model, elastics)


def change_constraint(
    this_constr: optlang.Constraint,
    new_expr: optlang.Constraint,
    model,
    sloppy: bool = False,
):
    """Change a constraint expression in the solver."""
    lb = this_constr.lb
    ub = this_constr.ub
    name = this_constr.name
    # Remove former constraint to override it
    model.remove(this_constr)
    new_cons = model.interface.Constraint(name=name, expression=new_expr, ub=ub, lb=lb)
    # Add the new variant
    model.add(new_cons, sloppy=sloppy)


def elastic_upper_relaxation(
    original_model: cobra.Model, elastic_candidates: List[str]
) -> Set:
    r"""Convert constrains to elastic constraints until the problem is feashible.

    It assumes that the elastic candidates are all subject to a <= constraint.

    Based on Brown and Graves, 1975:
    Q <- original problem
    E <- set of variables, candidates to be relaxes
    V <- set of variables that relax the problem
    relax(e_i) <- Add constraitn e - v
    Z <- \sum_{v_i \in V} v_i
    IIS <- minimal set of infeashible variables

    1. relax(e_i) for each e_i in E
    2. min Z, s.t. E, Q, v_i >= 0 for all V
    3. Get R <- v_i in V > 0
    4. If infeashible:
        STOP
    5. Else:
        5.1. E = E - R
        5.2. IIS = IIS & R
        go to 1
    6. return IIS
    """
    status = "optimal"
    current_candidates = set(elastic_candidates)
    iss = set()
    n_last_iss = -1
    while status == "optimal" and n_last_iss != len(iss):
        n_last_iss = len(iss)
        model = original_model.copy()
        # 5.1
        current_candidates -= iss
        objective_vars = []
        # 1.
        for e_id in current_candidates:
            e = model.constraints[e_id]
            v = model.problem.Variable(f"V_{e_id}", lb=0, ub=1000)
            new_expr = e.expression + v
            change_constraint(e, new_expr, model.solver)
            objective_vars.append(v)
        # 2.
        model.objective = model.problem.Objective(
            sympy.Add(*(objective_vars + [-model.objective.expression])),
            direction="min",
        )
        model.slim_optimize()
        status = model.solver.status
        # 5.2.
        iss |= {
            # 3.
            v.name[2:]
            for v in objective_vars
            if abs(v.primal) > 0 and v.name.startswith("V_")
        }
    return iss


def apply_upper_relaxation(
    original_model: cobra.Model, to_relax: List[str]
) -> Tuple[Model, Set]:
    """Relax the problem by relaxing `to_relax` upper (x <=) contraints."""
    model = original_model.copy()
    objective_vars = []
    for e_id in to_relax:
        e = model.constraints[e_id]
        v = model.problem.Variable(f"V_{e_id}", lb=0, ub=1000)
        new_expr = e.expression + v
        change_constraint(e, new_expr, model.solver)
        objective_vars.append(v)
    model.objective -= sympy.Add(*objective_vars)
    model.slim_optimize()
    return model, {v.name[2:] for v in objective_vars if abs(v.primal) > 0}


def top_shadow_prices(solution: cobra.Solution, top: int = 1) -> pd.DataFrame:
    """Rank them from most to least sensitive in the model reactions.

    Parameters
    ----------
    solution: cobra.Solution
        The usual Solution object returned by model.optimize().
    top: int
        The number of metabolites to be returned.

    Returns
    -------
    shadow_pr: pd.Series
        Top shadow prices, ranked.
    """
    shadow_pr = solution.shadow_prices
    return shadow_pr.sort_values()[:top]


def relax_proteomics_greedy(
    model: Model, minimal_growth: float
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

    # relax growth constraint
    minimal_growth *= 1.05

    # while the model cannot grow to the desired level, remove the protein with
    # the highest shadow price:
    prots_to_remove = []
    n_prots_with_concentration = len(
        [
            prot
            for prot in model.proteins
            if prot.concentration and not isnan(prot.concentration)
        ]
    )
    while new_growth_rate < minimal_growth and n_prots_with_concentration > 0:
        # get most influential protein in model:
        top_protein = top_shadow_prices(prot_solution)
        top_protein = model.proteins.get_by_id(top_protein.index[0])

        # update data: append protein to list, remove from current dataframe and
        # increase the corresponding upper bound to +1000:
        prots_to_remove.append(top_protein)
        print(top_protein.id)
        top_protein.add_concentration(None)
        top_protein.upper_bound = 1000
        n_prots_with_concentration -= 1

        # re-compute solution:
        solution, prot_solution = model.optimize()
        # if solution.objective_value == new_growth_rate:  # the algorithm is stuck
        #    break
        new_growth_rate = (
            solution.objective_value
            if solution.objective_value and not isnan(solution.objective_value)
            else 0
        )

    # update growth rate if optimization was not successful:
    if new_growth_rate < minimal_growth:
        print(
            f"Minimal growth was not reached! "
            f"Final growth of the model: {new_growth_rate}"
        )

    return new_growth_rate, prots_to_remove
