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

"""Protein manipulations for SBML compliance."""

import warnings
from collections import Counter
from enum import Enum
from typing import List, Optional, Union

from cobra.core.group import Group

from geckopy.model import Model, Protein


class EcStoichiometry(Enum):
    """Fashion of stoichiometric coeffients for enzymes (GECKO < or >= 3.0)."""

    KCAT = "KCAT"
    MW_KCAT = "MW_KCAT"


def group_proteins(
    model: Model, protein_list: Optional[List[Union[str, Protein]]] = None
):
    """Include proteins in `proteins` in the SBML group :code:`Protein`.

    It will create the :code:`Protein` SBML if it was not present already.

    Parameters
    ----------
    model: geckopy.Model
    proteins: Optional[List]
        list of proteins to include in the protein group (`str` or
        :code:`geckopy.Protein`). Default: include all proteins that are not
        already in the protein group.
    """
    if "Protein" not in model.groups:
        model.groups.add(Group("Protein", kind="classification"))
    proteins_to_add = (
        [prot.id for prot in model.proteins]
        if protein_list is None
        else [prot.id if isinstance(prot, Protein) else prot for prot in protein_list]
    )
    model.groups.get_by_id("Protein").add_members(
        [model.proteins.get_by_id(prot_id) for prot_id in proteins_to_add]
    )


def annotate_gene_protein_rules(model: Model):
    """Point gene to proteins and viceversa using gene annotations.

    EXPERIMENTAL: This is a helper to ease the gene to protein conversion but
    more thought should be put in how to make the GPR into something typed in the
    SBML doc that relates to the actual data structures in the model.
    """
    warnings.warn(
        "This function is experimental and subject to change in the future.",
        stacklevel=2,
    )
    # add empty attributes for the proteins for consistency
    for prot in model.proteins:
        prot.gene = None
    for gene in model.genes:
        gene.protein = None
        if "uniprot" in gene.annotation:
            prot_id = gene.annotation["uniprot"]
            prot = None
            if prot_id in model.proteins:
                prot = model.proteins.get_by_id(prot_id)
            elif f"prot_{prot_id}" in model.proteins:
                prot = model.proteins.get_by_id(f"prot_{prot_id}")
            if prot is not None:
                gene.protein = prot
                prot.gene = gene
    for prot in model.proteins:
        if prot.gene is None:
            # try to relate a gene to a protein using the reactions as a proxy
            gene_counter = Counter(
                (
                    gene
                    for reac in prot.reactions
                    # a reaction might have more than one genes and proteins
                    if len(reac.genes) == 1
                    for gene in reac.genes
                )
            )
            if len(gene_counter) > 0:
                gene = gene_counter.most_common(1)[0][0]
                gene.protein = prot
                prot.gene = gene
