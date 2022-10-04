from copy import copy, deepcopy

from cobra.core.dictlist import DictList
from multitfa.core import Thermo_met, thermo_reaction, tmodel

import geckopy


class ThermoProtReaction(thermo_reaction):
    """
    Class representation of thermo reaction Object. We calculate the required thermodynamic constraints for performing tMFA. To do the constraints, we need Gibbs energy of reaction and transport.

    Parameters
    ----------
    cobra_rxn : cobra.core.Reaction
        Cobra reaction object, to copy the attributes from. We copy metabolites and genes.
    updated_model : core.tmodel, optional
        tmodel object, with updated thermo properties, by default None

    """

    def __init__(
        self,
        cobra_rxn,
        updated_model=None,
    ):
        self._model = updated_model
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for attr, value in cobra_rxn.__dict__.items():
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = copy(value)

        self._metabolites = {}
        for met, stoic in cobra_rxn._metabolites.items():
            # this is the only change; also account for proteins
            new_met = (
                self.model.metabolites.get_by_id(met.id)
                if met in self.model.metabolites
                else self.model.proteins.get_by_id(met.id)
            )
            self._metabolites[new_met] = stoic
            new_met._reaction.add(self)
        self._genes = set()
        for gene in cobra_rxn._genes:
            new_gene = self.model.genes.get_by_id(gene.id)
            self._genes.add(new_gene)
            new_gene._reaction.add(self)


class ThermoProtModel(tmodel):
    def __init__(
        self,
        model,
        Exclude_list=[],
        tolerance_integral=1e-9,
        compartment_info=None,
        membrane_potential=None,
        exclude_metabolites=[],
    ):

        self.compartment_info = compartment_info
        self.membrane_potential = membrane_potential

        do_not_copy_by_ref = {
            "metabolites",
            "reactions",
            "proteins",
            "genes",
            "notes",
            "annotation",
        }
        for attr in model.__dict__:
            if attr not in do_not_copy_by_ref:
                self.__dict__[attr] = model.__dict__[attr]

        self.metabolites = DictList()
        do_not_copy_by_ref = {"_reaction", "_model"}
        for metabolite in model.metabolites:
            new_met = Thermo_met(
                metabolite=metabolite,
                updated_model=self,
            )
            self.metabolites.append(new_met)

        self.genes = DictList()

        for gene in model.genes:
            new_gene = gene.__class__(None)
            for attr, value in gene.__dict__.items():
                if attr not in do_not_copy_by_ref:
                    new_gene.__dict__[attr] = (
                        copy(value) if attr == "formula" else value
                    )
            new_gene._model = self
            self.genes.append(new_gene)
        self.proteins = DictList()
        for protein in model.proteins:
            new_prot = Thermo_met(
                metabolite=protein,
                updated_model=self,
            )
            # proteins do not participate in dGf calculations
            self.proteins.append(new_prot)

        self.reactions = DictList()
        do_not_copy_by_ref = {"_model", "_metabolites", "_genes"}
        for reaction in model.reactions:
            # this is custom to make the reaction aware of proteins
            new_reaction = ThermoProtReaction(
                cobra_rxn=reaction,
                updated_model=self,
            )
            self.reactions.append(new_reaction)

        try:
            self._solver = deepcopy(model.solver)
            # Cplex has an issue with deep copies
        except Exception:  # pragma: no cover
            self._solver = copy(model.solver)  # pragma: no cover

        self.Exclude_list = Exclude_list
        self.solver.configuration.tolerances.integrality = tolerance_integral
        self._var_update = False
