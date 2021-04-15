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

"""Set session level fixtures."""

from os.path import dirname, join
from typing import List

import cobra
import pytest

import geckopy


config = cobra.Configuration()


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


@pytest.fixture(scope="function")
def thermodb():
    """Path to thermodynamics database from pytfa."""
    return join(dirname(__file__), "data", "thermo_data.thermodb")


@pytest.fixture(scope="function")
def mnx():
    """Path to xref reduced xref file."""
    return join(dirname(__file__), "data", "chem_xref_seedset.tsv")


@pytest.fixture(scope="session")
def compartment_data():
    """Path to compartment data from pytfa."""
    return join(dirname(__file__), "data", "compartment_data.json")


@pytest.fixture(scope="session")
def fva_targets() -> List[str]:
    """Path to compartment data from pytfa."""
    return [
        "SUCASPtpp",
        "arm_NDPK1",
        "SUCCt1pp",
        "FUMt1pp",
        "SUCFUMtpp",
        "ASPtpp",
        "arm_ALATA_L",
        "GAPDNo1",
        "ENONo1",
        "GLCtex_copy1",
        "TPINo1",
        "PGINo1",
        "ICDHyrNo1",
        "MDHNo1",
        "G6PDH2rNo1",
        "GHMT2rNo1",
        "RPENo1",
        "MTHFDNo1",
        "MTHFCNo1",
        "CBMKrNo1",
        "EX_glc__D_e",
    ]


@pytest.fixture(scope="session")
def protein_list() -> List[str]:
    """Proteins to constrain as a pool set."""
    return [
        "prot_P45578",
        "prot_P45578",
        "prot_P0A796",
        "prot_P06999",
        "prot_P23721",
        "prot_P15254",
        "prot_P0A9B2",
    ]
