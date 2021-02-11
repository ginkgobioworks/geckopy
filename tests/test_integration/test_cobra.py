"""Ensure that ecgem is consistent with cobrapy."""

import pandas as pd

import ecgem


def test_unconstrained_ec_model_is_cobra_model(ec_model, cobra_model):
    """Check that unconstrained ec_model returns the same maximum as the plain model."""
    assert round(cobra_model.slim_optimize(), 4) == round(ec_model.slim_optimize(), 4)


def test_constrained_ec_model_is_not_cobra_model(cobra_model, experimental_copy_number):
    """Check that constrained ec_model returns different maximum than the plain model."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = ecgem.experimental.from_copy_number(
        cobra_model.copy(),
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    assert round(cobra_model.slim_optimize(), 4) != round(ec_model.slim_optimize(), 4)


def test_from_cobrapy_works(cobra_model):
    """Generate ec_gem from cobrapy_model."""
    ec_model = ecgem.Model(cobra_model)
    assert len(ec_model.proteins) == 1259
