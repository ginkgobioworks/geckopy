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

"""Test Input/Output capabilities."""
import os

import pytest

import geckopy
from geckopy.experimental.molecular_weights import extract_proteins


def test_read_geckopy_from_file(path_ecoli_core):
    """Read model directly from file."""
    model = geckopy.io.read_sbml_ec_model(path_ecoli_core)
    assert len(model.proteins) == 55


def test_copy_geckopy(ec_model_core):
    """Check that deepcopy works."""
    copied = ec_model_core.copy()
    assert len(copied.proteins) == len(ec_model_core.proteins)
    assert len(copied.reactions) == len(ec_model_core.reactions)
    assert len(copied.metabolites) == len(ec_model_core.metabolites)


def test_parsing_captures_naming_convention(dummy_ec_model):
    """Check proteins rely on the naming convention prot_UNIPROT are parsed."""
    assert dummy_ec_model.proteins.query("prot_P0A805")


def test_parsing_captures_protein_group(dummy_ec_model):
    """Check members of Protein group are parsed as proteins."""
    assert dummy_ec_model.groups.query("Protein")
    assert dummy_ec_model.proteins.query("prot_P0A825")
    assert dummy_ec_model.proteins.query("dummy_prot")


def test_protein_parsing_does_not_get_normal_metabolites(dummy_ec_model):
    """Check normal metabolites are not parsed as proteins."""
    assert not dummy_ec_model.proteins.query("normal_met")
    assert dummy_ec_model.metabolites.query("normal_met")
    mets = set(dummy_ec_model.metabolites)
    prots = set(dummy_ec_model.proteins)
    assert mets ^ prots == mets | prots


def test_serialized_model_grows(slim_solution_core, ec_model_core):
    """Check that deserialized model grows at the same rate."""
    geckopy.io.write_sbml_ec_model(ec_model_core, "_tmpfull.xml")
    redeserialized = geckopy.io.read_sbml_ec_model(
        "_tmpfull.xml", hardcoded_rev_reactions=False
    )
    assert pytest.approx(redeserialized.slim_optimize()) == pytest.approx(
        slim_solution_core
    )
    os.remove("_tmpfull.xml")


def test_writing_with_mw_kcat_produces_same_solution(slim_solution_core, ec_model_core):
    """Check that deserialized model with GECKO 3.0 grows at the same rate."""
    df = extract_proteins(ec_model_core)
    for row in df.itertuples():
        ec_model_core.proteins.get_by_id(row[2]).mw = row[3]
    geckopy.io.write_sbml_ec_model(
        ec_model_core,
        "_tmpfull_mwkcat.xml",
        ec_stoichiometry=geckopy.io.EcStoichiometry.MW_KCAT,
    )
    redeserialized = geckopy.io.read_sbml_ec_model(
        "_tmpfull_mwkcat.xml",
        hardcoded_rev_reactions=False,
        ec_stoichiometry=geckopy.io.EcStoichiometry.MW_KCAT,
    )
    assert pytest.approx(redeserialized.slim_optimize()) == pytest.approx(
        slim_solution_core
    )
    os.remove("_tmpfull_mwkcat.xml")


def test_serialized_model_has_concentrations(dummy_ec_model):
    """Check that concentrations are properly saved on SBML serialization."""
    dummy_ec_model.proteins.prot_P0A825.concentration = 123
    geckopy.io.write_sbml_ec_model(dummy_ec_model, "_tmp.xml")
    redeserialized = geckopy.io.read_sbml_ec_model(
        "_tmp.xml", hardcoded_rev_reactions=False
    )
    assert redeserialized.proteins.prot_P0A825.concentration == 123
    os.remove("_tmp.xml")


def test_proteins_are_grouped_on_write(dummy_ec_model):
    """Check that grouped proteins not following naming are properly handled."""
    dummy_ec_model.add_proteins([geckopy.Protein("my_unconventional_protein")])

    assert (
        dummy_ec_model.proteins.prot_P0A805
        not in dummy_ec_model.groups.get_by_id("Protein").members
    )
    geckopy.io.write_sbml_ec_model(
        dummy_ec_model, "_tmp_auto_grouping.xml", group_untyped_proteins=True  # default
    )
    # proteins that were not grouped but follow the conventions are not grouped
    assert (
        dummy_ec_model.proteins.prot_P0A805
        not in dummy_ec_model.groups.get_by_id("Protein").members
    )
    redeserialized = geckopy.io.read_sbml_ec_model(
        "_tmp_auto_grouping.xml", hardcoded_rev_reactions=False
    )
    assert (
        redeserialized.proteins.my_unconventional_protein.id
        == "my_unconventional_protein"
    )
    os.remove("_tmp_auto_grouping.xml")


def test_grouped_proteins_are_correctly_deserialized(ec_model_core):
    """Check that deepcopy works."""
    copied = ec_model_core.copy()
    copied.reactions.NH4t.add_protein("W", 60)
    geckopy.io.write_sbml_ec_model(
        copied, "_tmp_with_prot.xml", group_untyped_proteins=True  # default
    )
    model = geckopy.io.read_sbml_ec_model("_tmp_with_prot.xml")
    assert pytest.approx(model.proteins.W.kcats["NH4t"]) == 60
    os.remove("_tmp_with_prot.xml")


def test_gene_and_proteins_point_to_each_other(ec_model_core):
    """Check that annotating genes and proteins point to each other."""
    geckopy.io.standard.annotate_gene_protein_rules(ec_model_core)
    assert ec_model_core.genes.get_by_id(
        "b1241"
    ).protein == ec_model_core.proteins.get_by_id("prot_P0A9Q7")
    assert (
        ec_model_core.proteins.get_by_id("prot_P0A9Q7")
        == ec_model_core.genes.get_by_id("b1241").protein
    )
    # all proteins in the EC core model are in a gene
    assert sum(gene.protein is not None for gene in ec_model_core.genes) == len(
        ec_model_core.proteins
    )
