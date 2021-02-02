"""Ensure that ecgem is consistent with cobrapy."""

import ecgem
import pandas as pd


def test_unconstrained_ec_model_is_cobra_model(ec_model, cobra_model):
    """Check that unconstrained ec_model returns the same maximum as the plain model."""
    assert round(cobra_model.slim_optimize(), 4) == round(ec_model.slim_optimize(), 4)


def test_constrained_ec_model_is_not_cobra_model(
    ec_model, cobra_model, path_eciML1515, experimental_copy_number
):
    """Check that unconstrained ec_model returns the same maximum as the plain model."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = ecgem.experimental.from_copy_number(
        cobra_model,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    assert round(cobra_model.slim_optimize(), 4) != round(ec_model.slim_optimize(), 4)
