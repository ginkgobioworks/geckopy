"""Extend cobra.Reaction to point to proteins."""

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
