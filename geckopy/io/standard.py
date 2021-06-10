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

from typing import List, Optional, Union

from cobra.core.group import Group

from geckopy.model import Model, Protein


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
