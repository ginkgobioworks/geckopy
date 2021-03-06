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

"""API for proteins."""

import hashlib
import logging
import re
from math import isinf, isnan
from typing import Dict, FrozenSet, Tuple, Union
from warnings import warn

import optlang
from cobra import Configuration, Metabolite, Object, Reaction
from cobra.exceptions import OptimizationError
from cobra.util.context import resettable
from cobra.util.solver import (
    check_solver_status,
    linear_reaction_coefficients,
    set_objective,
)
from cobra.util.util import format_long_string


LOGGER = logging.getLogger(__name__)
UNIPROT_PATTERN = re.compile(
    r"(?:prot_)?[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
config = Configuration()


class Protein(Object):
    """Representation of an enzyme.

    A protein sets an upper bound to a set of reactions, given the kcat and
    concentration. Adapted from `cobra.Reaction`.

    In terms of the inner LP model, `Proteins` populates both variables (as
    pseudorreactions with the mentioned upper_bound) and constraints (as metabolites).

    Attributes
    ----------
    id: str
        identifier of protein, should be an Uniprot ID or "prot_<Uniprot ID>"
    name: str
        human readable name
    concentration: float
        parsed from `initialAmount` of `Species` in the SBML specification.
    kcats: Kcats
        1 / stoichimetry coefficients of its reactions.
    mw: float
        TODO: parsed from initalParameters/calculate from formula
    contribution: float
        value of flux of variable, only accessible after optimizing the model.
    lower_bound: float
        should be 0
    upper_bound: float
        concentration * kcat if there is a concentration, mmw if is part of the
        pool constraint (unmeasured proteins) or 1000.
    formula: str
    charge: float
    """

    def __init__(
        self,
        id: Union[str, Metabolite] = None,
        concentration: float = None,
        molecular_weight: float = 0.0,
        name: str = "",
    ):
        """Initialize with id."""
        self.concentration = concentration
        self.mw = molecular_weight
        self._flux = None
        self.lower_bound = 0
        self._reaction = set()
        self._ub = config.upper_bound
        self.compartment = "c"
        self.__metabolites = {}
        self._annotation = {}
        if isinstance(id, Metabolite):
            self.from_metabolite(id)
        else:
            Object.__init__(self, id, name)
            self.formula = ""
            self.charge = 0.0
        self.kcats = Kcats(self)
        self._reaction.add(self)

    def _set_id_with_model(self, value):
        """Overwrite parent id setter to change constraints and variables."""
        if value in self.model.proteins:
            raise ValueError(f"This model already constains a protein '{value}'")
        forward_variable = self.forward_variable
        reverse_variable = self.reverse_variable
        self.model.constraints[self.id].name = value
        self._id = value
        self.model.proteins._generate_index()
        forward_variable.name = self.id
        reverse_variable.name = self.reverse_id

    def from_metabolite(self, met: Metabolite):
        """Initialize `Protein` from `Metabolite`; i.e., when reading from SBML."""
        if not UNIPROT_PATTERN.match(met.id):
            LOGGER.warning(f"Metabolite {met.id} does not use Uniprot ID.")
        Object.__init__(self, met.id, met.name)
        self.compartment = met.compartment
        self._reaction = met._reaction
        self.charge = met.charge
        self.formula = ""
        self._annotation = met._annotation

    def update_variable_bounds(self):
        """Sync object bounds with inner model variable bounds."""
        if self.model is None:
            return
        # We know that `lb <= ub`.
        if self.lower_bound > 0:
            self.forward_variable.set_bounds(
                lb=None if isinf(self.lower_bound) else self.lower_bound,
                ub=None if isinf(self.upper_bound) else self.upper_bound,
            )
            self.reverse_variable.set_bounds(lb=0, ub=0)
        elif self.upper_bound < 0:
            self.forward_variable.set_bounds(lb=0, ub=0)
            self.reverse_variable.set_bounds(
                lb=None if isinf(self.upper_bound) else -self.upper_bound,
                ub=None if isinf(self.lower_bound) else -self.lower_bound,
            )
        else:
            self.forward_variable.set_bounds(
                lb=0, ub=None if isinf(self.upper_bound) else self.upper_bound
            )
            self.reverse_variable.set_bounds(
                lb=0, ub=None if isinf(self.lower_bound) else -self.lower_bound
            )

    def suscribe_to_pool(self, mmw: float):
        """Change internal pseudorreaction to have the common pool as reactant."""
        if not hasattr(self, "_model"):
            warn(
                f"Cannot set kcat of Protein {self.id} which does not"
                f" belong to any model."
            )
            return
        if not hasattr(self._model, "common_protein_pool"):
            warn(
                f"Cannot set kcat of Protein {self.id} whose Model does not"
                f" have a protein pool. Try `model.add_pool()` first."
            )
            return
        self._model.constraints[
            self._model.common_protein_pool.id
        ].set_linear_coefficients(
            {
                self.forward_variable: -mmw,
                self.reverse_variable: mmw,
            }
        )
        self._model.common_protein_pool._reaction.add(self)
        self.metabolites = {self._model.common_protein_pool: -mmw}

    def unsuscribe_to_pool(self):
        """Change internal pseudorreaction to have the common pool as reactant."""
        if hasattr(self, "_model"):
            if hasattr(self._model, "common_protein_pool"):
                self._model.constraints[
                    self._model.common_protein_pool.id
                ].set_linear_coefficients(
                    {
                        self.forward_variable: 0,
                        self.reverse_variable: 0,
                    }
                )
                self._model.common_protein_pool._reaction.remove(self)
                self.metabolites = {}

    @property
    def concentration(self) -> float:
        r"""Get upper bounds as [E] (conventionally in $\frac{mmol}{gDW}$).

        [E] multiplied by the kcat (expressed in the reaction stoichiometry as
        1/kcat) yields $\frac{mmol}/{gDW h}$.

        Taken from [Benjam??n J S??nchez et al., 2016]
        (https://www.embopress.org/doi/full/10.15252/msb.20167411).
        """
        return self._concentration

    @concentration.setter
    @resettable
    def concentration(self, value):
        self.add_concentration(value)

    @property
    def upper_bound(self) -> float:
        r"""Get upper bounds as [E] (conventionally in $\frac{mmol}{gDW}$).

        [E] multiplied by the kcat (expressed in the reaction stoichiometry as
        1/kcat) yields $\frac{mmol}/{gDW h}$.

        Taken from [Benjam??n J S??nchez et al., 2016]
        (https://www.embopress.org/doi/full/10.15252/msb.20167411).
        """
        return (
            self._ub
            if self.concentration is None
            or isinf(self.concentration)
            or isnan(self.concentration)
            else self.concentration
        )

    @upper_bound.setter
    @resettable
    def upper_bound(self, value: float):
        self._ub = value
        self.update_variable_bounds()

    def add_concentration(self, value: float):
        """Add concentration value.

        It will unsuscribe the protein to the common protein pool, if suscribed.
        """
        self._concentration = value
        self.unsuscribe_to_pool()
        self.update_variable_bounds()

    @property
    def reverse_id(self) -> str:
        """Generate the id of reverse_variable from the reaction's id."""
        return "_".join(
            (self.id, "reverse", hashlib.md5(self.id.encode("utf-8")).hexdigest()[0:5])
        )

    @property
    def contribution(self) -> str:
        """Get primal value (analogous to flux) in the most recent solution."""
        return self.flux

    @property
    def flux(self) -> float:
        """Get flux value in the most recent solution.

        Flux is the primal value of the corresponding variable in the model.

        Warnings
        --------
        * Accessing reaction fluxes through a `Solution` object is the safer,
          preferred, and only guaranteed to be correct way. You can see how to
          do so easily in the examples.
        * Reaction flux is retrieved from the currently defined
          `self._model.solver`. The solver status is checked but there are no
          guarantees that the current solver state is the one you are looking
          for.
        * If you modify the underlying model after an optimization, you will
          retrieve the old optimization values.

        Raises
        ------
        RuntimeError
            If the underlying model was never optimized beforehand or the
            reaction is not part of a model.
        OptimizationError
            If the solver status is anything other than 'optimal'.
        AssertionError
            If the flux value is not within the bounds.

        """
        try:
            check_solver_status(self._model.solver.status)
            return self.forward_variable.primal - self.reverse_variable.primal
        except AttributeError:
            raise RuntimeError(f"Protein '{self.id}' is not part of a model.")
        # Due to below all-catch, which sucks, need to reraise these.
        except (RuntimeError, OptimizationError) as err:
            raise OptimizationError(err)
        # Would love to catch CplexSolverError and GurobiError here.
        except Exception as err:
            raise OptimizationError(
                f"Likely no solution exists. Original solver message: {err}."
            )

    @property
    def forward_variable(self) -> optlang.Variable:
        """Get `optlang.Variable` representing the forward flux.

        Returns
        -------
        optlang.interface.Variable
            An optlang variable for the forward flux or None if reaction is
            not associated with a model.

        """
        if self.model is not None:
            return self.model.variables[self.id]
        else:
            return None

    @property
    def reverse_variable(self) -> optlang.Variable:
        """Get `optlang.Variable` representing the reverse flux.

        Returns
        -------
        optlang.interface.Variable
            An optlang variable for the reverse flux or None if reaction is
            not associated with a model.

        """
        if self.model is not None:
            return self.model.variables[self.reverse_id]
        else:
            return None

    @property
    def objective_coefficient(self) -> float:
        """Get the coefficient for this reaction in a linear objective (float).

        Assuming that the objective of the associated model is summation of
        fluxes from a set of reactions, the coefficient for each reaction
        can be obtained individually using this property. A more general way
        is to use the `model.objective` property directly.
        """
        return linear_reaction_coefficients(self.model, [self]).get(self, 0)

    @objective_coefficient.setter
    def objective_coefficient(self, value: float):
        if self.model is None:
            raise AttributeError("Cannot assign objective to a missing model.")
        if self.flux_expression is not None:
            set_objective(self.model, {self: value}, additive=True)

    @property
    def bounds(self) -> Tuple[float, float]:
        """Get or set the bounds directly from a tuple.

        Convenience method for setting upper and lower bounds in one line
        using a tuple of lower and upper bound. Invalid bounds will raise an
        AssertionError.
        When using a `HistoryManager` context, this attribute can be set
        temporarily, reversed when the exiting the context.
        """
        return self.lower_bound, self.upper_bound

    @bounds.setter
    @resettable
    def bounds(self, value: (float, float)):
        lower, upper = value
        # Validate bounds before setting them.
        Reaction._check_bounds(lower, upper)
        self.lower_bound = lower
        self.upper_bound = upper

    @property
    def _metabolites(self):
        """Get metabolites."""
        return {**self.__metabolites, self: 1}

    @_metabolites.setter
    def _metabolites(self, metabolites: Dict):
        """Get metabolites."""
        self.__metabolites = metabolites

    @property
    def metabolites(self) -> Dict:
        """Get metabolite of the protein pseudoreaction as the protein itself."""
        return self._metabolites

    @metabolites.setter
    def metabolites(self, metabolites: Dict):
        self._metabolites = metabolites

    @property
    def reactions(self) -> FrozenSet[Reaction]:
        """Retrieve immutable private reactions property."""
        return frozenset(reac for reac in self._reaction if isinstance(reac, Reaction))

    @property
    def model(self):
        """Retrieve the model the reaction is a part of."""
        return self._model if hasattr(self, "_model") else None

    @property
    def reactants(self):
        """Return a list of reactants for the reaction."""
        return [k for k, v in self._metabolites.items() if v < 0]

    @property
    def products(self):
        """Return a list of products for the reaction."""
        return [k for k, v in self._metabolites.items() if v >= 0]

    def __setstate__(self, state):
        """Set attribute to point to model, mimicking `cobra.Reaction`."""
        # These are necessary for old pickles which store attributes
        # which have since been superceded by properties.
        if "reaction" in state:
            state.pop("reaction")
        if "gene_reaction_rule" in state:
            state["_gene_reaction_rule"] = state.pop("gene_reaction_rule")
        if "lower_bound" in state:
            state["_lower_bound"] = state.pop("lower_bound")
        if "upper_bound" in state:
            state["_upper_bound"] = state.pop("upper_bound")

        self.__dict__.update(state)
        for x in state["_metabolites"]:
            x._model = self._model
            x._reaction.add(self)

    def __str__(self) -> str:
        """Print str representation as id."""
        return f"Protein {self.id}"

    def _repr_html_(self) -> str:
        return f"""
        <table>
            <tr>
                <td><strong>Protein identifier</strong></td><td>{self.id}</td>
            </tr><tr>
                <td><strong>Name</strong></td><td>{format_long_string(self.name, 200)}
                </td>
            </tr><tr>
                <td><strong>Memory address</strong></td>
                <td>0x0%x{id(self)}</td>
            </tr><tr>
            </tr><tr>
                <td><strong>Concentration</strong></td><td>{self.concentration}</td>
            </tr><tr>
                <td><strong>Upper bound</strong></td><td>{self.upper_bound}</td>
            </tr><tr>
            </tr><tr>
                <td><strong>Mw</strong></td><td>{self.mw}</td>
            </tr><tr>
                <td><strong>In {len(self.reactions)} reaction(s)</strong></td><td>
                    {format_long_string(", ".join(
                        f'{r.id} ({kcat:.2f})' for r, kcat in self.kcats.items()
                    ), 200)}
                </td>
            </tr>
        </table>
        """


class Kcats:
    """Interface to modify kcats of a protein.

    The user interacts with kcats in 1/s, translated to h to the model in the
    form stoichiometry coefficients of protein participating in a given reaction.

    Example
    -------
    .. doctest::

        # dictionary of Reaction to kcat in s
        >>> model.proteins.prot_P0A796.kcats
        # the user inputs the kcat in 1/s
        >>> model.proteins.prot_P0A796.kcats["PFKNo1"] = 1 / 30
        # the corresponding stoichiometry value is in -h
        >>> ec_model.reactions.PFKNo1.metabolites[model.proteins.prot_P0A796] == -1/120

    """

    def __init__(self, prot: Protein):
        """Initialize with kcats from the model."""
        self._protein = prot
        self._update()

    def _update(self):
        """Build the map on the fly."""
        self._reac_to_kcat = {
            reac: -1 / reac.metabolites[self._protein] / 3600
            for reac in self._protein.reactions
        }

    def __getitem__(self, key: Union[Reaction, str]) -> float:
        """Return the kcat of the protein in the Reaaction `key`."""
        self._update()
        if self._model_warn(key, "get"):
            return
        if isinstance(key, str):
            return self._reac_to_kcat[self._protein.model.reactions.get_by_id(key)]
        else:
            return self._reac_to_kcat[key]

    def __setitem__(self, key: Union[Reaction, str], val: float):
        """Assing kcat as `val` in the given `key` Reaction (as 1/kcat)."""
        # TODO: wrong
        if self._model_warn(key, "set"):
            return
        model = self._protein._model
        if isinstance(key, str):
            reac = model.reactions.get_by_id(key)
        else:
            reac = key
        coefficient = -1 / (3600 * val)
        reac._metabolites[self._protein] = coefficient
        # update corresponding linear coefficients in solver
        non_splitted_reaction = (
            list(
                model.constraints[self._protein.id]
                .get_linear_coefficients([reac.reverse_variable])
                .values()
            )[0]
            < 0
        )
        model.constraints[self._protein.id].set_linear_coefficients(
            {
                reac.forward_variable: coefficient,
                reac.reverse_variable: coefficient
                if non_splitted_reaction
                else -coefficient,
            }
        )

    def _model_warn(self, key: Union[Reaction, str], action: str):
        if not hasattr(self._protein, "model"):
            warn(
                f"Cannot set kcat of Protein {self._protein.id} which does not"
                f" belong to any model (reaction: {key})."
            )
            return True

    def __iter__(self):
        """Iterate inner dict."""
        self._update()
        return self._reac_to_kcat.__iter__()

    def items(self):
        """Iterate inner items."""
        self._update()
        return self._reac_to_kcat.items()

    def keys(self):
        """Return inner dict's keys."""
        self._update()
        return self._reac_to_kcat.keys()

    def values(self):
        """Return inner dict's values."""
        self._update()
        return self._reac_to_kcat.values()

    def __repr__(self):
        """Represent inner dict."""
        self._update()
        return self._reac_to_kcat.__repr__()
