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
from geckopy.flux_analysis import (
    flux_variability_analysis,
    get_protein_bottlenecks,
    get_protein_usage_by_reaction_rate,
    protein_variability_analysis,
    rate_kcat_concentration,
)


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


def test_fva(ec_model_core, fva_targets):
    """Test that fva returns the expected results."""
    df = flux_variability_analysis(ec_model_core)
    assert ((df.maximum - df.minimum) > 1e-3).sum() == 38


def test_pva(ec_model_core):
    """Test that fva returns the expected results."""
    df = protein_variability_analysis(ec_model_core)
    assert ((df.maximum - df.minimum) > 1e-3).sum() == 27


def test_usage_x_reaction_rate_is_consistent(ec_model_core):
    """Check that the protein usage of top proteins increases with Glc uptake."""
    fluxes = get_protein_usage_by_reaction_rate(
        ec_model_core, "EX_glc__D_e", [-0.1, -1.0, -5.0], "Glc uptake"
    )
    assert (
        fluxes.loc[fluxes["Glc uptake"] == -0.1, "fluxes"].sum()
        < fluxes.loc[fluxes["Glc uptake"] == -1.0, "fluxes"].sum()
        < fluxes.loc[fluxes["Glc uptake"] == -5.0, "fluxes"].sum()
    )


def test_expected_bottleneck(ec_model_core):
    """Check that the protein usage of top proteins increases with Glc uptake."""
    for prot in ec_model_core.proteins:
        prot.mw = 33000
    ec_model_core.constrain_pool(0.00448, 0.65, 1)
    bottlenecks = get_protein_bottlenecks(ec_model_core, 10)
    assert (bottlenecks.protein == "prot_P27306").any()


def test_kcat_concentration_rate(ec_model_core):
    """Check that the protein usage of top proteins increases with Glc uptake."""
    kcat_x_flux = rate_kcat_concentration(
        ec_model_core, "prot_P0A9P0", "PDH", [10, 100, 300, 639]
    )
    assert kcat_x_flux[0][1] < kcat_x_flux[1][1] < kcat_x_flux[2][1]
