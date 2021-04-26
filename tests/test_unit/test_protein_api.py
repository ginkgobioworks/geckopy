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


def test_all_proteins_are_parsed(ec_model_core):
    """Test .proteins interface."""
    assert len(ec_model_core.proteins) == 55


def test_protein_point_to_reactions(ec_model_core):
    """Test protein.reactions interface."""
    assert (
        ec_model_core.reactions.get_by_id("PDH")
        in ec_model_core.proteins.prot_P0A9P0.reactions
    )
    assert len(ec_model_core.proteins.prot_P0A9P0.reactions) == 2


def test_reaction_point_to_proteins(ec_model_core):
    """Test protein.reactions.REAC_ID.proteins interface."""
    assert (
        ec_model_core.proteins.get_by_id("prot_P0A9P0")
        in ec_model_core.reactions.PDH.proteins
    )


def test_protein_pseudoreactions_are_not_model_reactions(ec_model):
    """Ensure that reactions of legacy EC model are not prot pseudorreactions."""
    assert not ec_model.reactions.query("prot_P75905_exchange")


def test_protein_are_not_model_metabolites(ec_model_core):
    """Ensure that model.metabolites is not polluted with proteins."""
    assert not ec_model_core.metabolites.query("prot_P0A9P0")


def test_kcats_retrieve_right_coefficients(ec_model_core):
    """Test kcat interface."""
    a_prot = ec_model_core.proteins.get_by_id("prot_P0A9P0")
    assert int(a_prot.kcats["PDH"]) == int(37.90)
    a_prot.kcats["PDH"] = 2
    assert ec_model_core.reactions.get_by_id("PDH").metabolites[a_prot] == -1 / (
        2 * 3600
    )
