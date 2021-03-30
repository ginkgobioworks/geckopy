"""Model class that extends `cobra.Model` to account for enzyme constraints."""
import logging
from collections import defaultdict
from copy import copy, deepcopy
from functools import partial
from math import isnan
from typing import Dict, Iterator, Optional, Tuple, Union

import cobra
import pandas as pd
from cobra import Metabolite
from cobra.core.dictlist import DictList
from cobra.util.context import get_context
from optlang.symbolics import Zero

from .protein import Protein
from .reaction import Reaction


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
        # TODO: separate proteins from cobra.Reactions
        The last obtained solution from optimizing the model.

    """

    def __setstate__(self, state: Dict):
        """Make sure all cobra.Objects in the model point to the model."""
        self.__dict__.update(state)
        for y in ["reactions", "genes", "metabolites", "proteins"]:
            for x in getattr(self, y):
                x._model = self
        if not hasattr(self, "name"):
            self.name = None

    def __init__(
        self,
        id_or_model: Union[str, cobra.Model] = None,
        name: str = None,
        hardcoded_rev_reactions: bool = True,
    ):
        """Initialize model."""
        # TODO: refactor from cobra.Model so that internal matrix does not get
        # repopulated 100 times
        self.proteins = DictList()
        self.hardcoded_rev_reactions = hardcoded_rev_reactions
        if isinstance(id_or_model, cobra.Model) and not hasattr(
            id_or_model, "proteins"
        ):
            super().__init__(id_or_model.id, name)
            self.from_cobra(id_or_model)
        else:
            super().__init__(id_or_model, name)

    def from_cobra(self, model):
        """Create from cobra model."""
        g_proteins = [met for met in model.metabolites if met.id.startswith("prot_")]
        # groups
        if model.groups.query("Protein"):
            g_proteins += model.groups.Protein.members
        g_proteins = set(g_proteins)
        g_prot_id = {prot.id for prot in g_proteins}
        self._solver = model.solver
        self.add_proteins([Protein(prot) for prot in g_proteins])
        self.add_metabolites(
            [met for met in model.metabolites if met.id not in g_prot_id]
        )
        self.add_reactions(
            [reac for reac in model.reactions if not reac.id.startswith("prot_")]
        )
        self.__setstate__(self.__dict__)
        self._populate_solver(self.reactions, self.metabolites, self.proteins)

    def get_total_measured_proteins(self) -> float:
        """Sum of all `Proteins` in the model that has a concentration."""
        measured = sum(prot.concentration for prot in self.proteins)
        return 0 if isnan(measured) or not measured else measured

    def constrain_pool(
        self, p_total, sigma_saturation_factor, fn_mass_fraction_unmeasured_matched
    ):
        """Constrain the draw reactions for the unmeasured (common protein pool) proteins.

        Adapted from [geckopy]
        (https://github.com/SysBioChalmers/GECKO/blob/master/geckopy/geckopy/gecko.py#L184)

        Proteins without their own protein pool are collectively constrained by the
        common protein pool. Remove protein pools for all proteins that don't have
        measurements, along with corresponding draw reactions, and add these to
        the common protein pool and reaction.

        Parameters
        ----------
        p_total: float
            measured total protein fraction in cell in g protein / gDW
        sigma_saturation_factor: float
            part of proteome that can be used by metabolism
        fn_mass_fraction_unmeasured_matched: float
            TODO: add convenience function to handle this
            sum of the product of average abundances of unmesured proteins
            (from, e.g., paxDB) times their molecular weight.
        """
        unmeasured_prots = [
            prot
            for prot in self.proteins
            if prot.concentration and ~isnan(prot.concentration)
        ]
        if not unmeasured_prots:
            return
        self.add_pool()
        # if there are no proteins measured, the calculations are simplified
        fs_matched_adjusted = p_total * sigma_saturation_factor
        if len(unmeasured_prots) != len(self.proteins):
            p_measured = self.get_total_measured_proteins()
            # * section 2.5.1
            # 1. and 2. introduce `prot_pool` and exchange reaction done in __init__
            # 3. limiting total usage with the unmeasured amount of protein
            f_mass_fraction_measured_matched_to_total = (
                fn_mass_fraction_unmeasured_matched / (1 - p_measured / p_total)
            )
            fs_matched_adjusted = (
                (p_total - p_measured)
                * f_mass_fraction_measured_matched_to_total
                * sigma_saturation_factor
            )
        self.protein_pool_exchange.bounds = 0, fs_matched_adjusted

        m_weigths = [prot.mw for prot in self.proteins if prot.mw]
        average_mmw = sum(m_weigths) / len(m_weigths) / 1000.0
        for protein in unmeasured_prots:
            mmw = protein.mw / 1000 if protein.mw else average_mmw
            protein.suscribe_to_pool(mmw)

    def add_pool(self, reac_id: str = "prot_pool_exchange", met_id: str = "prot_pool"):
        """Add the reaction and metabolite to constraint the common pool of proteins.

        Parameters
        ----------
        read_id: str
            id of the common protein pool pseudorreaction to add.
        met_id: str
            id of the common protein pool metabolite to add.

        """
        try:
            self.common_protein_pool = self.metabolites.get_by_id(met_id)
        except KeyError:
            self.common_protein_pool = Metabolite(met_id)
        try:
            self.protein_pool_exchange = self.reactions.get_by_id(reac_id)
        except KeyError:
            self.protein_pool_exchange = Reaction(reac_id)
            self.protein_pool_exchange.add_metabolites({self.common_protein_pool: 1.0})
            self.add_reactions([self.protein_pool_exchange])

    def copy(self):
        """Provide a partial 'deepcopy' of the Model.

        All of the Metabolite, Gene, and Reaction objects are created anew but
        in a faster fashion than deepcopy.
        Enzyme constrained changes: also deepcopy proteins.
        """
        new = self.__class__()
        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "proteins",
            "genes",
            "notes",
            "annotation",
            "groups",
        }
        for attr in self.__dict__:
            if attr not in do_not_copy_by_ref:
                new.__dict__[attr] = self.__dict__[attr]
        new.notes = deepcopy(self.notes)
        new.annotation = deepcopy(self.annotation)

        new.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in self.metabolites:
            new_met = metabolite.__class__()
            for attr, value in metabolite.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_met.__dict__[attr] = copy(value) if attr == "formula" else value
            new_met._model = new
            new.metabolites.append(new_met)

        new.proteins = DictList()
        for protein in self.proteins:
            new_prot = protein.__class__()
            for attr, value in protein.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_prot.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_prot._model = new
            new.proteins.append(new_prot)

        new.genes = DictList()
        for gene in self.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = new
            new.genes.append(new_gene)

        new.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in self.reactions:
            new_reaction = reaction.__class__()
            for attr, value in reaction.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_reaction.__dict__[attr] = copy(value)
            new_reaction._model = new
            new.reactions.append(new_reaction)
            # update awareness
            for metabolite, stoic in reaction._metabolites.items():
                if metabolite in new.proteins:
                    new_met = new.proteins.get_by_id(metabolite.id)
                    new_reaction._metabolites[new_met] = stoic
                    new_met._reaction.add(new_reaction)
                else:
                    # regular met
                    new_met = new.metabolites.get_by_id(metabolite.id)
                    new_reaction._metabolites[new_met] = stoic
                    new_met._reaction.add(new_reaction)
            for gene in reaction._genes:
                new_gene = new.genes.get_by_id(gene.id)
                new_reaction._genes.add(new_gene)
                new_gene._reaction.add(new_reaction)

        new.groups = DictList()
        do_not_copy_by_ref = {"_model", "_members"}
        # Groups can be members of other groups. We initialize them first and
        # then update their members.
        for group in self.groups:
            new_group = group.__class__(group.id)
            for attr, value in group.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_group.__dict__[attr] = copy(value)
            new_group._model = new
            new.groups.append(new_group)
        for group in self.groups:
            new_group = new.groups.get_by_id(group.id)
            # update awareness, as in the reaction copies
            new_objects = []
            for member in group.members:
                if isinstance(member, Metabolite):
                    new_object = new.metabolites.get_by_id(member.id)
                elif isinstance(member, Reaction):
                    new_object = new.reactions.get_by_id(member.id)
                elif isinstance(member, cobra.Gene):
                    new_object = new.genes.get_by_id(member.id)
                elif isinstance(member, Protein):
                    new_object = new.proteins.get_by_id(member.id)
                elif isinstance(member, cobra.core.Group):
                    new_object = new.genes.get_by_id(member.id)
                else:
                    raise TypeError(
                        "The group member {!r} is unexpectedly not a "
                        "metabolite, reaction, gene, nor another "
                        "group.".format(member)
                    )
                new_objects.append(new_object)
            new_group.add_members(new_objects)

        try:
            new._solver = deepcopy(self.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            new._solver = copy(self.solver)  # pragma: no cover

        # it doesn't make sense to retain the context of a copied model so
        # assign a new empty context
        new._contexts = list()

        return new

    def _populate_solver(
        self,
        reaction_list: Iterator[Reaction],
        metabolite_list: Iterator[Metabolite] = None,
        protein_list: Iterator[Protein] = None,
    ):
        """Populate attached solver with LP problem given reactions + proteins.

        Note that proteins are added both as Constraints and Variables.
        """
        constraint_terms = defaultdict(lambda: defaultdict(float))
        to_add = []
        metabolite_list = metabolite_list if metabolite_list is not None else []
        protein_list = protein_list if protein_list is not None else []
        if metabolite_list or protein_list:
            for met in metabolite_list + protein_list:
                if met.id not in self.constraints:
                    # this condition only applies when passing `protein_list`
                    to_add += [self.problem.Constraint(Zero, name=met.id, lb=0, ub=0)]
        self.add_cons_vars(to_add)

        for reaction in reaction_list + protein_list:
            reaction_id = reaction.id

            if reaction_id not in self.variables:
                forward_variable = self.problem.Variable(reaction_id)
                reverse_variable = self.problem.Variable(reaction.reverse_id)
                self.add_cons_vars([forward_variable, reverse_variable])
            else:
                try:
                    reaction = self.reactions.get_by_id(reaction_id)
                except KeyError:
                    reaction = self.proteins.get_by_id(reaction_id)
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
                # if a reaction is reversible, the protein is used on both
                # directions
                rev = (
                    reaction.reversibility
                    # it could be a protein
                    if hasattr(reaction, "reversibility")
                    else False
                )
                constraint_terms[constraint][reverse_variable] = (
                    coeff
                    if not self.hardcoded_rev_reactions
                    and metabolite in self.proteins
                    and rev
                    else -coeff
                )

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
        self._populate_solver([], None, pruned)

    def add_reactions(self, reaction_list: Iterator[Reaction]):
        """Add reactions to the model.

        Reactions with identifiers identical to a reaction already in the
        model are ignored.
        The change is reverted upon exit when using the model as a context.
        Enzyme Constrained changes: avoid adding proteins as metabolites.

        Parameters
        ----------
        reaction_list : list
            A list of `cobra.Reaction` objects
        """

        def existing_filter(rxn):
            if rxn.id in self.reactions:
                LOGGER.warning(f"Ignoring reaction '{rxn.id}' since it already exists.")
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
                if metabolite not in target:
                    self.add_metabolites(metabolite)

                # A copy of the metabolite exists in the model, the reaction
                # needs to point to the metabolite in the model.
                else:
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

    def optimize(
        self, objective_sense: Optional[str] = None, raise_error: bool = False
    ) -> Tuple[cobra.Solution, cobra.Solution]:
        """Optimize the model using flux balance analysis.

        Parameters
        ----------
        objective_sense : {None, 'maximize' 'minimize'}, optional
            Whether fluxes should be maximized or minimized. In case of None,
            the previous direction is used.
        raise_error : bool
            If true, raise an OptimizationError if solver status is not
             optimal.

        Notes
        -----
        Only the most commonly used parameters are presented here.  Additional
        parameters for cobra.solvers may be available and specified with the
        appropriate keyword argument.

        """
        solution = super().optimize(objective_sense, raise_error)
        solution_prot = cobra.core.get_solution(
            self, reactions=self.proteins, metabolites=self.proteins
        )
        solution_prot.contributions = solution_prot.fluxes
        # workaround for wrong shadow prices of proteins
        try:
            prot_index = [prot.id for prot in self.proteins]
            solution_prot.shadow_prices = pd.Series(
                index=prot_index,
                data=[self.constraints[prot.id].dual for prot in prot_index],
                name="shadow_prices",
            )
        except Exception:
            pass

        return solution, solution_prot

    def add_boundary(
        self,
        metabolite: Union[cobra.Metabolite, Protein],
        type: str = "exchange",
        reaction_id: Optional[str] = None,
        lb: Optional[float] = None,
        ub: Optional[float] = None,
        sbo_term: Optional[str] = None,
    ) -> Reaction:
        """Add a boundary reaction for a given metabolite.

        Enzyme constraint changes: return an geckopy.Reaction.
        """
        rxn = super().__init__(metabolite, type, reaction_id, lb, ub, sbo_term)
        return Reaction(rxn)
