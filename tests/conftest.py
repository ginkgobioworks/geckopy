"""Set session level fixtures."""

from os.path import dirname, join

import cobra
import pytest


config = cobra.Configuration()
config.solver = "glpk"


@pytest.fixture(scope="session")
def eciML1515():
    """Load from cobrapy."""
    model = cobra.io.read_sbml_model(
        join(dirname(__file__), "data", "eciML1515.xml.gz")
    )
    return model


@pytest.fixture(scope="session")
def experimental_copy_number():
    return join(dirname(__file__), "data", "ecoli_proteomics_schmidt2016S5.tsv")


@pytest.fixture(scope="session")
def experimental_mmol_gDW():
    return join(dirname(__file__), "data", "mmol_gdW_protemics.csv")
