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

"""COBRA flux Methods to work with/around EC model nuisances."""

import logging
import re
from itertools import chain
from math import isnan
from multiprocessing import Pool
from typing import List, Optional, Tuple

import cobra
import pandas as pd
from cobra.util.solver import set_objective
from numpy import zeros
from optlang.symbolics import Zero
from tqdm import tqdm

from .model import Model


LOGGER = logging.getLogger(__name__)
config = cobra.Configuration()
REV_PATTERN = re.compile(r"(No\d)|$")
UNREV_PATTERN = re.compile(r"(.+)_REV(No\d)?")
__all__ = [
    "flux_variability_analysis",
    "protein_variability_analysis",
    "pfba_protein",
    "get_protein_bottlenecks",
    "get_protein_usage_by_reaction_rate",
    "rate_kcat_concentration",
]


class UserInputModelClash(Exception):
    """Error caused by inconsistencies between the model and the user input."""

    pass


def _apply_rev(reac_id: str, model: cobra.Model) -> str:
    rev_id = re.sub(REV_PATTERN, r"_REV\1", reac_id)
    return rev_id if rev_id in model.reactions else None


def _reverse_rev(reac_id: str, model: cobra.Model) -> str:
    rev_id = re.sub(UNREV_PATTERN, r"\1\2", reac_id)
    return rev_id if rev_id != reac_id and rev_id in model.reactions else None


def _reverse_arm(reac_id: str, model: cobra.Model) -> str:
    rev_id = "arm_" + re.sub(REV_PATTERN, r"", reac_id)
    return rev_id if rev_id != reac_id and rev_id in model.reactions else None


def _get_objective_reaction(model: cobra.Model) -> cobra.Reaction:
    return next(var for var in model.objective.variables if var.name in model.reactions)


def _sanitize_input_list(fixed_reactions: Optional[List[str]] = None) -> List:
    return (
        []
        if fixed_reactions is None
        else [fixed_reactions]
        if isinstance(fixed_reactions, str)
        else fixed_reactions
    )


def _init_worker(model, target="reactions"):
    """Initialize a global model object for multiprocessing."""
    global _model
    global _target
    _model = model
    _target = target


def isnan_none(x: float) -> bool:
    return isnan(x) or x is None


def _fva_step(reaction_id):
    """Run a maximization and a minization.

    Modified from cobrapy to account for possibly duplicated (`_REV`) reactions
    in the EC model.
    """
    global _model
    global _target
    target = _model.__getattribute__(_target)
    reac = target.get_by_id(reaction_id)
    rev_id = _apply_rev(reaction_id, _model)
    rev_reac = target.get_by_id(rev_id) if rev_id else None
    set_objective(_model, {reac: 1})
    if rev_reac is not None:
        # if reac with duplicate is optimized, block counterpart
        prev_bounds = rev_reac.bounds
        rev_reac.bounds = 0, 0
        ub = _model.slim_optimize()
        rev_reac.bounds = prev_bounds
        # and use the reverse for the minimum (blocking forward)
        prev_bounds = reac.bounds
        set_objective(_model, {rev_reac: 1}, False)
        reac.bounds = 0, 0
        lb = -_model.slim_optimize()
        reac.bounds = prev_bounds
    else:
        ub = _model.slim_optimize()
        set_objective(_model, {reac: -1}, False)
        lb = _model.slim_optimize()
    # handle infeasible case
    if isnan_none(lb) or isnan_none(ub):
        lb = 0.0 if isnan_none(lb) else lb
        ub = 0.0 if isnan_none(ub) else ub
        LOGGER.warning(
            "Could not get flux for reaction %s,  setting "
            "it to 0. This is usually due to numerical instability.",
            reaction_id,
        )
    return reaction_id, lb, ub


def fix_reaction_to_min(model: cobra.Model, reac: cobra.Reaction) -> Tuple[str, float]:
    """Fix a reaction to its minimum."""
    rev_id = _apply_rev(reac.id, model)
    rev_reac = model.reactions.get_by_id(rev_id) if rev_id else None
    with model:
        if rev_reac is not None:
            prev_bounds = reac.bounds
            set_objective(model, {rev_reac: 1})
            reac.bounds = 0, 0
            lb = model.slim_optimize()
            reac.bounds = prev_bounds
            reac_to_fix = rev_id
        else:
            set_objective(model, {reac: -1})
            lb = model.slim_optimize()
            reac_to_fix = reac.id
    return reac_to_fix, lb


def flux_variability_analysis(
    ec_model: Model,
    fixed_reactions: Optional[List[str]] = None,
    ignored_reactions: Optional[List[str]] = None,
    n_proc: Optional[int] = config.processes,
    inplace: bool = False,
) -> pd.DataFrame:
    r"""Flux variability analysis for EC models.

    This method is detailed in the Supplementary material of [Domenzain et al.,
    2021](https://www.biorxiv.org/content/10.1101/2021.03.05.433259v1.full.pdf).

    For a fair comparison of flux distributions, EC models require some adjustments:

    #. Since there are duplicated reactions (for each direction),
       the complementary reaction must be blocked to remove artificial flux variation.
    #. Some reactions might be fixed (the glucose exchange reaction in the
       paper) to compare them with non-EC models.
    #. The final reported flux is
       :math:`v_{max}, v_{max}^{rev} \forall v \in \text{reactions}`

    Additionally, just the combined reaction `arm_` is reported for isozymes.

    See Also
    --------
    :py:func:`cobrapy:cobra.flux_analysis.flux_variability_analysis`


    Parameters
    ----------
    model: geckopy.Model
    fixed_reactions: Optional[List[str]]
        List of reactions to be fixed (as in the second point in the description)
    ignored_reactions: Optional[List[str]]
        List of reactions to be ignored
    n_proc: int
        Number of processes to use. Default: the number of logical CPUs.
    inplace: bool
        Whether to copy the model to then apply the changes

    Returns
    -------
    fva_result: pandas.DataFrame
        A data frame with reaction identifiers as the index and two columns:
        - maximum: indicating the highest possible flux
        - minimum: indicating the lowest possible flux

    """
    fixed_reactions = _sanitize_input_list(fixed_reactions)
    ignored_reactions = _sanitize_input_list(ignored_reactions)

    model = ec_model if inplace else ec_model.copy()
    fixed_objective = model.slim_optimize()
    _get_objective_reaction(model).lb = fixed_objective

    reac_ids = [
        reac.id
        for reac in model.reactions
        # REV reactions will be taken into account in its forward counterpart
        if _reverse_rev(reac.id, model) is None and
        # keep arm reactions (gathering isozymes)
        _reverse_arm(reac.id, model) is None and reac.id not in ignored_reactions
    ]
    fva_result = pd.DataFrame(
        {
            "minimum": zeros(len(reac_ids), dtype=float),
            "maximum": zeros(len(reac_ids), dtype=float),
        },
        index=reac_ids,
    )

    reac_ids = [reac for reac in reac_ids if reac not in fixed_reactions]
    for reac in fixed_reactions:
        if reac not in model.reactions:
            raise UserInputModelClash(f"{reac} to be fixed is not in model {model.id}.")
        reac_to_fix, lb = fix_reaction_to_min(model, model.reactions.get_by_id(reac))
        model.reactions.get_by_id(reac_to_fix).bounds = lb, lb
        fva_result.loc[reac, :] = [-lb, -lb] if reac != reac_to_fix else [lb, lb]
    LOGGER.debug(
        f"FVA started for a model with {len(model.reactions)} and"
        f"{len(reac_ids)} effective reactions"
    )

    chunk_size = len(reac_ids) // n_proc
    if n_proc > 1:
        with Pool(
            n_proc,
            initializer=_init_worker,
            initargs=(model,),
        ) as pool:
            for rxn_id, lb, ub in tqdm(
                pool.imap_unordered(_fva_step, reac_ids, chunksize=chunk_size),
                total=len(reac_ids),
                desc="FVA",
            ):
                fva_result.loc[rxn_id, :] = [lb, ub]
    else:
        _init_worker(model)
        for rxn_id, lb, ub in tqdm(
            map(_fva_step, reac_ids),
            total=len(reac_ids),
            desc="FVA",
        ):
            fva_result.loc[rxn_id, :] = [lb, ub]

    return fva_result


def annotate_protein_genes(sr: pd.Series, model: Model) -> pd.Series:
    """Build a series of genes mapping to a list of proteins."""
    return sr.apply(
        lambda x: ",".join(
            [
                ",".join([gene.id for gene in reac.genes])
                for reac in model.proteins.get_by_id(x).reactions
            ]
        )
    )


def get_protein_bottlenecks(model: Model, top: int = 10) -> pd.DataFrame:
    """Return `top` protein bottlenecks based on the shadow prices."""
    _, prots = model.optimize()
    df = prots.to_frame().sort_values("reduced_costs", ascending=False).head(top)
    df = df.reset_index().rename({"index": "protein"}, axis=1)
    df["gene"] = annotate_protein_genes(df.protein, model)
    return df


def protein_variability_analysis(
    ec_model: Model,
    protein_list: Optional[List[str]] = None,
    n_proc: Optional[int] = config.processes,
) -> pd.DataFrame:
    """Return the top used proteins for each reaction rate.

    Parameters
    ----------
    model: geckopy.Model
    protein_list: list[float]
    n_proc: Optional[int]

    Return
    ------
    pandas.DataFrame
        with columns 'minimum', 'maximum' and index as `Protein.id`

    """
    model = ec_model.copy()
    fix_solution = model.slim_optimize()
    _get_objective_reaction(model).bounds = fix_solution, fix_solution
    proteins = (
        protein_list[:]
        if protein_list is not None
        else [prot.id for prot in model.proteins]
    )
    fva_result = pd.DataFrame(
        {
            "minimum": zeros(len(proteins), dtype=float),
            "maximum": zeros(len(proteins), dtype=float),
        },
        index=proteins,
    )
    chunk_size = len(proteins) // n_proc
    if n_proc > 1:
        with Pool(
            n_proc,
            initializer=_init_worker,
            initargs=(model, "proteins"),
        ) as pool:
            for rxn_id, lb, ub in tqdm(
                pool.imap_unordered(_fva_step, proteins, chunksize=chunk_size),
                total=len(proteins),
                desc="FVA",
            ):
                fva_result.loc[rxn_id, :] = lb, ub
    else:
        _init_worker(model, "proteins")
        for rxn_id, lb, ub in tqdm(
            map(_fva_step, proteins),
            total=len(proteins),
            desc="FVA",
        ):
            fva_result.loc[rxn_id, :] = lb, ub
    return fva_result


def pfba_protein(
    model: Model, objective: Optional[str] = None, fraction_of_optimum: float = 1.0
) -> cobra.core.Solution:
    """Add pFBA objective.

    Add objective to minimize the summed flux of all proteins to the
    current objective.

    See Also
    -------
    :py:func:`cobrapy:cobra.flux_analysis.pfba`

    Parameters
    ----------
    model : cobra.Model
        The model to add the objective to
    objective :
        An objective to set in combination with the pFBA objective.
    fraction_of_optimum : float
        Fraction of optimum which must be maintained. The original objective
        reaction is constrained to be greater than maximal_value *
        fraction_of_optimum.

    Result
    ------
    solution: cobra.core.Solution

    """
    with model as m:
        if objective is not None:
            model.objective = objective
        if model.solver.objective.name == "_pfba_objective":
            raise ValueError("The model already has a pFBA objective.")
        cobra.util.fix_objective_as_constraint(model, fraction=fraction_of_optimum)
        reaction_variables = (
            (rxn.forward_variable, rxn.reverse_variable) for rxn in model.proteins
        )
        variables = chain(*reaction_variables)
        model.objective = model.problem.Objective(
            Zero, direction="min", sloppy=True, name="_pfba_objective"
        )
        model.objective.set_linear_coefficients({v: 1.0 for v in variables})
        m.slim_optimize(error_value=None)
        solution = cobra.core.solution.get_solution(m, reactions=model.proteins)
    return solution


def get_protein_usage_by_reaction_rate(
    model: Model,
    reaction: str,
    fix_fluxes: List[float],
    reaction_plot_name: str = "Reaction rate",
    top: int = 10,
) -> pd.DataFrame:
    """Return the top used proteins for each reaction rate.

    Parameters
    ----------
    model: geckopy.Model
    reaction: str
        id of the reaction whose flux will be fixed to each `fix_fluxes`.
    fix_fluxes: list[float]
        the model will be optimized by sequentally fixing the reaction bounds
        to these fluxes.
    top: int
        number of top proteins to return for each condition

    Return
    ------
    pandas.DataFrame
        with colums 'protein', 'fluxes', 'reduced_costs', `reaction_plot_name`, 'gene'

    Example
    -------
    We can make a bar plot with a slider for each uptake rate using
    `plotly express <https://plotly.com/python/plotly-express/>`_:

    .. doctest::

        >>> import plotly.express as px
        >>> model.constrain_pool(0.448, 0.65, 0.8)
        >>> fluxes = geckopy.flux_analysis.get_protein_usage_by_reaction_rate(
                model, "EX_glc__D_e", [-0.1, -0.2, -0.5, -1.0, -5.0, -10.0], "Glc uptake"
            )
        >>> # return percentages of the total protein pool
        >>> fluxes.fluxes *= 100 / 0.448
        >>> px.bar(
                fluxes, x="protein", y="fluxes", hover_data=["gene"], color="fluxes",
                animation_frame="Glc uptake",
                color_continuous_scale="TealGrn", template="simple_white"
            )
    """
    results = []
    prev_bounds = model.reactions.get_by_id(reaction).bounds
    for flux in fix_fluxes:
        model.reactions.get_by_id(reaction).bounds = (flux, flux)
        prot_usage = pfba_protein(model).to_frame()
        prot_usage[reaction_plot_name] = flux
        prot_usage = prot_usage.reset_index().rename({"index": "protein"}, axis=1)
        results.append(prot_usage)
    model.reactions.get_by_id(reaction).bounds = prev_bounds
    proteins = pd.concat(results).reset_index(drop=True)
    to_plot = (
        proteins.groupby(reaction_plot_name)
        .apply(lambda x: x.sort_values("fluxes", ascending=False).head(top))
        .reset_index(drop=True)
    )
    to_plot["gene"] = annotate_protein_genes(to_plot.protein, model)
    return to_plot.sort_values([reaction_plot_name, "fluxes"], ascending=False)


def rate_kcat_concentration(
    model: Model, protein: str, reaction: str, kcat_range: List[float]
):
    """Calculate the flux value for each kcat in `kcat_range`."""
    ec_model = model.copy()
    prot = ec_model.proteins.get_by_id(protein)

    def optimize_with_kcat(kcat):
        prot.kcats[reaction] = kcat
        return ec_model.slim_optimize()

    return [(optimize_with_kcat(kcat), kcat) for kcat in kcat_range]


def count_protein_saturation(model: Model):
    """Compute the number of proteins that are saturated."""
    return sum(
        1 for prot in model.proteins if abs(prot.upper_bound - prot.contribution) < 1e-7
    )


def get_average_protein_saturation(model: Model):
    """Compute the average upper bounds usage by the proteins in the model."""
    non_zero_prots = [prot for prot in model.proteins if prot.upper_bound]
    return sum(prot.contribution / prot.upper_bound for prot in non_zero_prots) / len(
        non_zero_prots
    )
