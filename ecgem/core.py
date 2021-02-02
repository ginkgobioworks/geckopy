"""Model class that extends `cobra.Model` to account for enzyme constraints."""
import logging
from collections import defaultdict
from functools import partial
from typing import Union, Iterator

import cobra
from cobra import Reaction, Metabolite
from cobra.core.dictlist import DictList
from cobra.util.context import get_context
from optlang.symbolics import Zero

from .protein import Protein

LOGGER = logging.getLogger(__name__)


class Model(cobra.Model):
    """Extension of cobra.Model providing an API for proteins in EC models.

    Attributes
    ----------
    reactions : DictList
        A DictList where the key is the reaction identifier and the value a Reaction.
    metabolites : DictList
        A DictList where the key is the metabolite identifier and the value a
        Metabolite.
    proteins : DictList
        A DictList where the key is the metabolite identifier and the value a Protein.
    genes : DictList
        A DictList where the key is the gene identifier and the value a Gene.
    groups : DictList
        A DictList where the key is the group identifier and the value a Group.
    solution : Solution
        The last obtained solution from optimizing the model.

    """

    def __setstate__(self, state):
        """Make sure all cobra.Objects in the model point to the model."""
        self.__dict__.update(state)
        for y in ["reactions", "genes", "metabolites", "proteins"]:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __init__(self, id_or_model: Union[str, cobra.Model] = None, name: str = None):
        """Initialize model."""
        super().__init__(id_or_model, name)
        if isinstance(id_or_model, str):
            self.proteins = DictList()

    def copy(self):
        """Copy (deep) from `cobra.Model` + proteins."""
        new = super().copy()
        # TODO: copy proteins
        new.proteins = DictList()
        return new

    def _populate_solver(
        self,
        reaction_list: Iterator[Reaction],
        metabolite_list: Iterator[Metabolite] = None,
        protein_list: Iterator[Protein] = None,
    ):
        """Populate attached solver with LP problem given reactions + proteins."""
        constraint_terms = defaultdict(lambda: defaultdict(float))
        to_add = []
        metabolite_list = metabolite_list if metabolite_list is not None else []
        protein_list = protein_list if protein_list is not None else []
        if metabolite_list or protein_list:
            for met in metabolite_list + protein_list:
                if met.id not in self.constraints:
                    to_add += [self.problem.Constraint(Zero, name=met.id, lb=0, ub=0)]
        self.add_cons_vars(to_add)

        for reaction in reaction_list + protein_list:
            # reaction_id = (
            #     reaction.id if reaction in reaction_list else f"{reaction.id}_ex"
            # )
            reaction_id = reaction.id

            if reaction_id not in self.variables:
                forward_variable = self.problem.Variable(reaction_id)
                reverse_variable = self.problem.Variable(reaction.reverse_id)
                self.add_cons_vars([forward_variable, reverse_variable])
            else:
                reaction = self.reactions.get_by_id(reaction_id)
                forward_variable = reaction.forward_variable
                reverse_variable = reaction.reverse_variable
            for metabolite, coeff in reaction.metabolites.items():
                if metabolite.id in self.constraints:
                    constraint = self.constraints[metabolite.id]
                else:
                    constraint = self.problem.Constraint(
                        Zero, name=metabolite.id, lb=0, ub=0
                    )
                    self.add_cons_vars(constraint, sloppy=True)
                constraint_terms[constraint][forward_variable] = coeff
                constraint_terms[constraint][reverse_variable] = -coeff

        self.solver.update()
        for reaction in reaction_list:
            self.reactions.get_by_id(reaction.id).update_variable_bounds()
        for protein in protein_list:
            self.proteins.get_by_id(protein.id).update_variable_bounds()
        for constraint, terms in constraint_terms.items():
            constraint.set_linear_coefficients(terms)

    def add_proteins(self, protein_list: Iterator[Protein]):
        """Add proteins to the model, in the same fashion as `.add_metabollites`."""
        if not hasattr(protein_list, "__iter__"):
            protein_list = [protein_list]
        if len(protein_list) == 0:
            return None
        # Take left difference
        pruned = [x for x in protein_list if x.id not in self.proteins]
        for prot in pruned:
            prot._model = self

        to_add = []
        for prot in pruned:
            if prot.id not in self.constraints:
                constraint = self.problem.Constraint(Zero, name=prot.id, lb=0, ub=0)
                to_add += [constraint]

        self.add_cons_vars(to_add)

        context = get_context(self)
        if context:
            context(partial(self.proteins.__isub__, pruned))
            for x in pruned:
                context(partial(setattr, x, "_model", None))

        self.proteins += pruned
        self._populate_solver([], None, self.proteins)

    def add_reactions(self, reaction_list):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.
        The change is reverted upon exit when using the model as a context.
        Enzyme Constrained: this have to be changed to avoid adding proteins
        as metabolites.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
        """

        def existing_filter(rxn):
            if rxn.id in self.reactions:
                LOGGER.warning(
                    "Ignoring reaction '%s' since it already exists.", rxn.id
                )
                return False
            return True

        # First check whether the reactions exist in the model.
        pruned = DictList(filter(existing_filter, reaction_list))

        context = get_context(self)

        # Add reactions. Also take care of genes and metabolites in the loop.
        for reaction in pruned:
            reaction._model = self
            # Build a `list()` because the dict will be modified in the loop.
            for metabolite in list(reaction.metabolites):
                target = (
                    self.proteins if metabolite in self.proteins else self.metabolites
                )
                # TODO: Should we add a copy of the metabolite instead?
                if metabolite not in target:
                    self.add_metabolites(metabolite)

                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
                    # FIXME: Modifying 'private' attributes is horrible.
                    stoichiometry = reaction._metabolites.pop(metabolite)

                    model_metabolite = target.get_by_id(metabolite.id)
                    reaction._metabolites[model_metabolite] = stoichiometry
                    model_metabolite._reaction.add(reaction)
                    if context:
                        context(partial(model_metabolite._reaction.remove, reaction))

            for gene in list(reaction._genes):
                # If the gene is not in the model, add it
                if not self.genes.has_id(gene.id):
                    self.genes += [gene]
                    gene._model = self

                    if context:
                        # Remove the gene later
                        context(partial(self.genes.__isub__, [gene]))
                        context(partial(setattr, gene, "_model", None))

                # Otherwise, make the gene point to the one in the model
                else:
                    model_gene = self.genes.get_by_id(gene.id)
                    if model_gene is not gene:
                        reaction._dissociate_gene(gene)
                        reaction._associate_gene(model_gene)

        self.reactions += pruned

        if context:
            context(partial(self.reactions.__isub__, pruned))

        # from cameo ...
        self._populate_solver(pruned)
