"""Tests for experimental data loading."""

import pandas as pd

from ecgem.experimental import from_copy_number, from_mmol_gDW


def test_model_from_copy_number_can_grow(cobra_model, experimental_copy_number):
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        cobra_model,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    assert ec_model.slim_optimize() > 0


def test_model_from_mmol_gDW_can_grow(cobra_model, experimental_mmol_gDW):
    processed_proteomics = pd.read_csv(experimental_mmol_gDW)
    ec_model = from_mmol_gDW(cobra_model, processed_proteomics)
    assert ec_model.slim_optimize() > 0
