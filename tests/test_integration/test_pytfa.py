"""Test main functionality."""
import os

import pytest
import pytfa

from pytfa.io import load_thermoDB

from geckopy.integration.pytfa import (
    prepare_gecko_to_thermo,
    add_dummy_protein_info,
    write_thermodb,
)


def test_integrated_model_works(ec_model, thermodb, mnx):
    """Test building block."""
    thermo_loaded = load_thermoDB(thermodb)
    prepare_gecko_to_thermo(ec_model, thermo_loaded, mnx)
    tmodel = pytfa.ThermoModel(thermo_loaded, ec_model)
    assert tmodel.slim_optimize() > 1e-2


def test_thermo_constrain_solution(ec_model, thermodb, mnx):
    """Check thermo model returns different solution than normal model."""
    sol = ec_model.optimize()[0]
    summed_sol = sol.fluxes.sum()
    thermo_loaded = load_thermoDB(thermodb)
    prepare_gecko_to_thermo(ec_model, thermo_loaded, mnx)
    thermo_model = pytfa.ThermoModel(thermo_loaded, ec_model)
    tsol = thermo_model.optimize()
    tsummed_sol = tsol.fluxes.sum()
    assert pytest.approx(tsummed_sol) != pytest.approx(summed_sol)


def test_thermo_with_protein_constrain(ec_model, thermodb, mnx):
    """Check ec-thermo model returns different solution than thermo model."""
    thermo_loaded = load_thermoDB(thermodb)
    prepare_gecko_to_thermo(ec_model, thermo_loaded, mnx)
    thermo_model = pytfa.ThermoModel(thermo_loaded, ec_model)
    tsol = thermo_model.slim_optimize()
    ec_model.proteins.prot_P0AEN1.add_concentration(2e-4)
    thermo_model = pytfa.ThermoModel(thermo_loaded, ec_model)
    tsol_ec_constrained = thermo_model.slim_optimize()
    assert pytest.approx(tsol) != pytest.approx(tsol_ec_constrained)


def test_write_thermodb(thermodb, dummy_ec_model):
    """Check that the thermodb can be written to a file."""
    thermodb = load_thermoDB(thermodb)
    add_dummy_protein_info(dummy_ec_model, thermodb)
    write_thermodb(thermodb, "thermodb.thermodb")
    changed_tdb = load_thermoDB("thermodb.thermodb")
    assert changed_tdb["metabolites"]["protein"]["deltaGf_std"] == 0
    # clean artifact
    os.remove("thermodb.thermodb")
