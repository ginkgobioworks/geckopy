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

"""Tests related to `.proteins` API."""

from math import isnan

from pytest import approx

import geckopy


def test_constraint_pool_changes_objective_value(slim_solution_core, ec_model_core):
    """Test .proteins interface."""
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model_core.proteins:
        prot.mw = 330
    ec_model_core.constrain_pool(2e-4, 0.8, 1.0)
    pool_solution = ec_model_core.slim_optimize()
    assert approx(slim_solution_core) != pool_solution
    assert approx(pool_solution) != 0.0 and not isnan(pool_solution)


def test_constraint_pool_set_changes_solution(
    slim_solution_core, ec_model_core, protein_list
):
    """Test .proteins interface."""
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model_core.proteins:
        prot.mw = 330
    ec_model_core.constrain_pool(0.0002, 0.065, 1.0, protein_list=protein_list)
    pool_solution = ec_model_core.slim_optimize()
    assert approx(slim_solution_core) != pool_solution
    assert approx(pool_solution) != 0.0 and not isnan(pool_solution)
    assert (
        len(ec_model_core.constraints["prot_pool"].variables)
        == len(protein_list) * 2 + 2
    )


def test_added_protein_modifies_solution(ec_model_core, slim_solution_core):
    """Test that adding a protein constrains the solution."""
    ec_model_core.reactions.CYTBD.add_protein(
        id="prot_INVNTD", kcat=0.3, concentration=2e-5
    )
    assert approx(ec_model_core.slim_optimize()) != slim_solution_core


def test_added_reaction_gathers_proteins(ec_model_core):
    """Add a reaction with a protein and check that it is correctly structured."""
    model = geckopy.Model("one_reaction_model")
    rxn = ec_model_core.reactions.GLUSy
    model.add_reactions([rxn])
    assert len(model.reactions) == 1
    assert len(model.metabolites) == 6
    assert {prot.id for prot in model.proteins} == {"prot_P09832", "prot_P09831"}
