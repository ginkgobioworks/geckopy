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


def test_constraint_pool_changes_objective_value(slim_solution, ec_model):
    """Test .proteins interface."""
    # mw are not included in the model, let's say all of them are 330 Da
    for prot in ec_model.proteins:
        prot.mw = 330
    ec_model.constrain_pool(2e-3, 0.8, 1.0)
    pool_solution = ec_model.slim_optimize()
    assert approx(slim_solution) != approx(pool_solution)
    assert approx(pool_solution) != 0.0 and not isnan(pool_solution)


def test_added_protein_modifies_solution(ec_model, slim_solution):
    """Test that adding a protein constrains the solution."""
    ec_model.reactions.CYTBO3_4ppNo1.add_protein(
        id="prot_INVNTD", kcat=0.3, concentration=2e-5
    )
    assert approx(ec_model.slim_optimize()) != approx(slim_solution)
