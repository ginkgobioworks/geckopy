"""API for proteins."""

import hashlib
import logging
import re
from math import isinf, isnan
from typing import Union

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
    kcat: float
        parsed from annotation of `Species`.
    mw: float
        TODO: parsed from initalParameters/calculate from formula
    flux: float
        value of flux of variable, only accessible after optimizing the model.
    lower_bound: float
        should be 0
    upper_bound: float
        concentration * kcat, it cannot be set directly
    formula: str
    charge: float
    """

    def __init__(
        self,
        id: Union[str, Metabolite] = None,
        concentration: float = None,
        kcat: float = 1.0,
        molecular_weight: float = 0.0,
    ):
        """Initialize with id."""
        self.concentration = concentration
        self.kcat = kcat
        self.mw = molecular_weight
        self._flux = None
        self.lower_bound = 0
        self._reaction = set()
        if isinstance(id, Metabolite):
            self.from_metabolite(id)
        else:
            self.id = id
            self.formula = ""
            self.charge = 0.0

    def from_metabolite(self, met: Metabolite):
        """Initialize `Protein` from `Metabolite`; i.e., when reading from SBML."""
        if not UNIPROT_PATTERN.match(met.id):
            LOGGER.warning(f"Metabolite {met.id} does not use Uniprot ID.")
        self.id = met.id
        self.name = met.name
        self._reaction = met._reaction
        self.charge = met.charge
        self.formula = met.formula

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

    @property
    def reverse_id(self):
        """Generate the id of reverse_variable from the reaction's id."""
        return "_".join(
            (self.id, "reverse", hashlib.md5(self.id.encode("utf-8")).hexdigest()[0:5])
        )

    @property
    def upper_bound(self):
        r"""Get upper bounds as [E] (conventionally in $\frac{mmol}{gDW}$).

        [E] multiplied by the kcat (expressed in the reaction stoichiometry as
        1/kcat) yields $\frac{mmol}/{gDW h}$.

        Taken from [Benjamín J Sánchez et al., 2016]
        (https://www.embopress.org/doi/full/10.15252/msb.20167411).
        """
        return (
            config.upper_bound
            if self.concentration is None
            or isinf(self.concentration)
            or isnan(self.concentration)
            else self.concentration
        )

    @upper_bound.setter
    @resettable
    def upper_bound(self, value):
        self.concentration = value
        self.update_variable_bounds()

    @property
    def flux(self):
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
    def forward_variable(self):
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
    def reverse_variable(self):
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
    def objective_coefficient(self):
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
    def bounds(self):
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
    def metabolites(self):
        """Get metabolite of the protein pseudoreaction as the protein itself."""
        return {self: 1}

    @property
    def reactions(self):
        """Retrieve immutable private reactions property."""
        return frozenset(self._reaction)

    @property
    def model(self):
        """Retrieve the model the reaction is a part of."""
        return self._model

    def __str__(self):
        """Print str representation as id."""
        return f"Protein {self.id}"

    def _repr_html_(self):
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
                <td><strong>In {len(self.reactions)} reaction(s)</strong></td><td>
                    {format_long_string(", ".join(r.id for r in self.reactions), 200)}
                </td>
            </tr>
        </table>
        """
