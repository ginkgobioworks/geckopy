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
from math import isnan
from multiprocessing import Pool
from typing import List, Optional, Tuple

import cobra
import pandas as pd
from cobra.util.solver import set_objective
from numpy import zeros
from tqdm import tqdm

from .model import Model


LOGGER = logging.getLogger(__name__)
config = cobra.Configuration()
REV_PATTERN = re.compile(r"(No\d)|$")
UNREV_PATTERN = re.compile(r"(.+)_REV(No\d)?")
__all__ = [
    "flux_variability_analysis",
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


def _init_worker(model):
    """Initialize a global model object for multiprocessing."""
    global _model
    _model = model


def isnan_none(x: float) -> bool:
    return isnan(x) or x is None


def _fva_step(reaction_id):
    """Run a maximization and a minization.

    Modified from cobrapy to account for possibly duplicated (`_REV`) reactions
    in the EC model.
    """
    global _model
    reac = _model.reactions.get_by_id(reaction_id)
    rev_id = _apply_rev(reaction_id, _model)
    rev_reac = _model.reactions.get_by_id(rev_id) if rev_id else None
    set_objective(_model, {reac: 1})
    if rev_reac is not None:
        # if reac with duplicate is optimized, block counterpart
        prev_bounds = rev_reac.bounds
        rev_reac.bounds = 0, 0
        ub = _model.slim_optimize()
        rev_reac.bounds = prev_bounds
        # and use the reverse for the minimum (blocking forward)
        prev_bounds = reac.bounds
        set_objective(_model, {rev_reac: 1})
        reac.bounds = 0, 0
        lb = -_model.slim_optimize()
        reac.bounds = prev_bounds
    else:
        ub = _model.slim_optimize()
        set_objective(_model, {reac: -1})
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

    Parametes
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
        fva_result.at[reac, :] = (-lb, -lb) if reac != reac_to_fix else (lb, lb)
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
                fva_result.at[rxn_id, :] = lb, ub
    else:
        _init_worker(model)
        for rxn_id, lb, ub in tqdm(
            map(_fva_step, reac_ids),
            total=len(reac_ids),
            desc="FVA",
        ):
            fva_result.at[rxn_id, :] = lb, ub

    return fva_result
