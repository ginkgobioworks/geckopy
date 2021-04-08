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

"""Test specialized flux analysis methods.

This tests are for GLPK and they are not consistent with CPLEX and gurobi.
However, CPLEX and gurobi are both consistent in the same solution.
"""
from geckopy import flux_variability_analysis


def test_fva_with_fixed_reactions(ec_model, fva_targets):
    """Test that fva with a fixed reaction returns the expected results."""
    df = flux_variability_analysis(
        ec_model,
        fixed_reactions="EX_glc__D_e",
        ignored_reactions=[
            reac.id for reac in ec_model.reactions if reac.id not in fva_targets
        ],
    )
    assert ((df.maximum - df.minimum) > 1e-3).sum() == 20


def test_fva(ec_model, fva_targets):
    """Test that fva returns the expected results."""
    df = flux_variability_analysis(
        ec_model,
        ignored_reactions=[
            reac.id for reac in ec_model.reactions if reac.id not in fva_targets
        ],
    )
    assert ((df.maximum - df.minimum) > 1e-3).sum() == 21
