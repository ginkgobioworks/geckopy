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


def test_constraint_pool_changes_objective_value(slim_solution, ec_model):
    """Test .proteins interface."""
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model.proteins:
        prot.mw = 330
    ec_model.constrain_pool(2e-3, 0.8, 1.0)
    pool_solution = ec_model.slim_optimize()
    assert approx(slim_solution) != approx(pool_solution)
    assert approx(pool_solution) != 0.0 and not isnan(pool_solution)


def test_constraint_pool_set_changes_solution(slim_solution, ec_model, protein_list):
    """Test .proteins interface."""
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model.proteins:
        prot.mw = 330
    ec_model.constrain_pool(0.0002, 0.65, 1.0, protein_list=protein_list)
    pool_solution = ec_model.slim_optimize()
    assert approx(slim_solution) != approx(pool_solution)
    assert approx(pool_solution) != 0.0 and not isnan(pool_solution)
    assert len(ec_model.constraints["prot_pool"].variables) == len(protein_list) * 2


def test_added_protein_modifies_solution(ec_model, slim_solution):
    """Test that adding a protein constrains the solution."""
    ec_model.reactions.CYTBO3_4ppNo1.add_protein(
        id="prot_INVNTD", kcat=0.3, concentration=2e-5
    )
    assert approx(ec_model.slim_optimize()) != approx(slim_solution)


def test_added_reaction_gathers_proteins(ec_model):
    """Add a reaction with a protein and check that it is correctly structured."""
    model = geckopy.Model("one_reaction_model")
    rxn = ec_model.reactions.PUACGAMtrNo1
    model.add_reaction(rxn)
    assert len(model.reactions) == 1
    assert len(model.metabolites) == 2
    assert {prot.id for prot in model.proteins} == {"prot_P75905", "prot_P69432"}
