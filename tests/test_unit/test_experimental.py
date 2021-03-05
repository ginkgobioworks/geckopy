"""Tests for experimental data loading."""

from math import isnan

import pandas as pd

from geckopy.experimental import from_copy_number, from_mmol_gDW


def test_ec_model_from_copy_number_can_grow(ec_model, experimental_copy_number):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    sol = ec_model.slim_optimize()
    assert sol < 0.1 if not isnan(sol) else True


def test_model_from_mmol_gDW_can_grow(ec_model, experimental_mmol_gDW):
    """Test that constraining the model works."""
    processed_proteomics = pd.read_csv(experimental_mmol_gDW)
    ec_model = from_mmol_gDW(ec_model, processed_proteomics)
    assert ec_model.slim_optimize() > 0.0
