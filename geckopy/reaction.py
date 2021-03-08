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

    def add_protein(self, id: str, kcat: float, concentration: Optional[float] = None):
        """Add a protein to the reaction.

        Parameters
        ==========
        id: str
        kcat: float (1/s)
        concentration: float
        """
        prot = Protein(id, kcat)
        coefficient = -1 / (kcat * 3600)
        prot.add_concentration(concentration)
        _id_to_metabolites = dict([(x.id, x) for x in self._metabolites])

        prot_id = str(prot.id)
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

            for metabolite, coefficient in self._metabolites.items():
                model.constraints[metabolite.id].set_linear_coefficients(
                    {
                        self.forward_variable: coefficient,
                        self.reverse_variable: -coefficient,
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
