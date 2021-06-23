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

"""Extend cobra.Reaction to point to proteins."""

from typing import Optional

import cobra

from .protein import Protein


class Reaction(cobra.Reaction):
    """geckopy.Reaction has an additional attribute `.proteins`."""

    @property
    def proteins(self):
        """Retrieve metabolites participating in the reaction that are proteins."""
        return {
            met: stoich
            for met, stoich in self._metabolites.copy().items()
            if isinstance(met, Protein)
        }

    def add_protein(
        self,
        id: str,
        kcat: float,
        concentration: Optional[float] = None,
        name: Optional[str] = None,
    ):
        """Add a protein to the reaction.

        Parameters
        ==========
        id: str
        kcat: float (1/s)
        concentration: float
        name: str
        """
        prot_id = str(id)
        name = name if name is not None else prot_id
        prot = Protein(prot_id, name=name)
        coefficient = -1 / (kcat * 3600)
        prot.add_concentration(concentration)
        _id_to_metabolites = {x.id: x for x in self._metabolites}

        if prot_id in _id_to_metabolites:
            self._metabolites[_id_to_metabolites[prot_id]] = coefficient
        else:
            # If the reaction is in a model, ensure we aren't using
            # a duplicate metabolite.
            if self._model:
                try:
                    prot = self._model.proteins.get_by_id(prot_id)
                except KeyError:
                    pass
            self._metabolites[prot] = coefficient
            # make the metabolite aware that it is involved in this
            # reaction
            prot._reaction.add(self)

        # from cameo ...
        model = self.model
        if model is not None:
            model.add_proteins([prot])

            # the protein is added to both sides of the reaction in case the
            # latter is reversible (it has to be consumed in both cases)
            model.constraints[prot.id].set_linear_coefficients(
                {
                    self.forward_variable: coefficient,
                    self.reverse_variable: coefficient,
                }
            )
        # context = get_context(self)
        # if context and reversibly:
        #     if combine:
        #         # Just subtract the metabolites that were added
        #         context(
        #             partial(
        #                 self.subtract_metabolites,
        #                 metabolites_to_add,
        #                 combine=True,
        #                 reversibly=False,
        #             )
        #         )
        #     else:
        #         # Reset them with add_metabolites
        #         mets_to_reset = {
        #             key: old_coefficients[model.metabolites.get_by_any(key)[0]]
        #             for key in iterkeys(metabolites_to_add)
        #         }

        #         context(
        #             partial(
        #                 self.add_metabolites,
        #                 mets_to_reset,
        #                 combine=False,
        #                 reversibly=False,
        #             )
        #         )

    def remove_protein(self, protein: Protein):
        """Remove protein from the reaction."""
        model = self.model
        if model is not None:
            model.constraints[protein.id].set_linear_coefficients(
                {self.forward_variable: 0, self.reverse_variable: 0}
            )
        protein._reaction.remove(self)
        self._metabolites.pop(protein)
