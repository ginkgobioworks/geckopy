import pytest
from optlang import available_solvers

from geckopy.integration.multitfa import ThermoProtModel


@pytest.mark.skipif(not (available_solvers["GUROBI"] or available_solvers["CPLEX"]))
def test_protein_constrain_affects_multitfa_solution(ec_model_core):
    """Check thermo model returns different solution when protein is constrained."""
    thermo_model = ThermoProtModel(ec_model_core.copy())
    tsol = thermo_model.slim_optimize()
    ec_model_core.proteins.prot_P25516.add_concentration(2e-5)
    thermo_model = ThermoProtModel(ec_model_core)
    tsol_prot = thermo_model.slim_optimize()
    assert pytest.approx(tsol) != tsol_prot
