"""Set session level fixtures."""

from os.path import dirname, join

import cobra
import pytest

import geckopy


config = cobra.Configuration()
config.solver = "glpk"


@pytest.fixture(scope="session")
def path_eciML1515():
    """Store path to model."""
    return join(dirname(__file__), "data", "eciML1515.xml.gz")


@pytest.fixture(scope="function")
def ec_model(path_eciML1515):
    """Load from cobrapy."""
    return geckopy.io.read_sbml_ec_model(path_eciML1515)


@pytest.fixture(scope="function")
def cobra_model(path_eciML1515):
    """Load from cobrapy."""
    return cobra.io.read_sbml_model(path_eciML1515)


@pytest.fixture(scope="session")
def dummy_ec_model():
    """Load from cobrapy."""
    return geckopy.io.read_sbml_ec_model(
        join(dirname(__file__), "data", "all_prot_encodings.xml")
    )


@pytest.fixture(scope="session")
def slim_solution(path_eciML1515):
    """Provide the objective value of the enzyme constrained model."""
    return geckopy.io.read_sbml_ec_model(path_eciML1515).slim_optimize()


@pytest.fixture(scope="session")
def experimental_copy_number():
    """Load exp data.

    TODO: add ref.
    """
    return join(dirname(__file__), "data", "ecoli_proteomics_schmidt2016S5.tsv")


@pytest.fixture(scope="session")
def experimental_mmol_gDW():
    """Load exp processed data.

    TODO: add ref.
    """
    return join(dirname(__file__), "data", "mmol_gdW_protemics.csv")
