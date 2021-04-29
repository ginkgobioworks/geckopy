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

"""Test main functionality."""
import os

import pytest
import pytfa
from pytfa.io import load_thermoDB

from geckopy.integration.pytfa import (
    adapt_gecko_to_thermo,
    get_thermo_coverage,
    translate_model_mnx_to_seed,
    write_thermodb,
)


def test_integrated_model_works(ec_model_core, thermodb, mnx, compartment_data):
    """Test building block."""
    thermodb = load_thermoDB(thermodb)
    compartment_data = pytfa.io.read_compartment_data(compartment_data)
    translate_model_mnx_to_seed(ec_model_core, thermodb, mnx)
    tmodel = adapt_gecko_to_thermo(ec_model_core, thermodb, compartment_data)
    assert get_thermo_coverage(tmodel) == 19


def test_thermo_constrain_solution(ec_model_core, thermodb, compartment_data, mnx):
    """Check thermo model returns different solution that normal model."""
    sol = ec_model_core.optimize()[0]
    summed_sol = sol.fluxes.sum()
    thermodb = load_thermoDB(thermodb)
    compartment_data = pytfa.io.read_compartment_data(compartment_data)
    translate_model_mnx_to_seed(ec_model_core, thermodb, mnx)
    tmodel = adapt_gecko_to_thermo(ec_model_core, thermodb, compartment_data)
    tsol = tmodel.optimize()
    tsummed_sol = tsol.fluxes.sum()
    assert pytest.approx(tsummed_sol) != summed_sol


def test_thermo_with_protein_constrain(ec_model_core, thermodb, compartment_data, mnx):
    """Check thermo model returns different solution that normal model."""
    thermodb = load_thermoDB(thermodb)
    compartment_data = pytfa.io.read_compartment_data(compartment_data)
    translate_model_mnx_to_seed(ec_model_core, thermodb, mnx)
    tmodel = adapt_gecko_to_thermo(ec_model_core, thermodb, compartment_data)
    tsol = tmodel.slim_optimize()
    ec_model_core.proteins.prot_P25516.add_concentration(2e-4)
    tmodel = adapt_gecko_to_thermo(ec_model_core, thermodb, compartment_data)
    tsol_ec_constrained = tmodel.slim_optimize()
    assert pytest.approx(tsol) != tsol_ec_constrained


def test_write_thermodb(thermodb):
    """Check that the thermodb can be written to a file."""
    thermodb = load_thermoDB(thermodb)
    thermodb["metabolites"]["protein"] = {
        "pKa": [7],
        "deltaGf_err": 0,
        "mass_std": 333.0,
        "struct_cures": {},
        "id": "protein",
        "nH_std": 12,
        "name": "protein",
        "formula": "",
        "deltaGf_std": 0,
        "error": "Nil",
        "charge_std": 0,
        "struct_cues": {},
    }
    write_thermodb(thermodb, "thermodb.thermodb")
    changed_tdb = load_thermoDB("thermodb.thermodb")
    # clean artifact
    os.remove("thermodb.thermodb")
    assert changed_tdb["metabolites"]["protein"]["deltaGf_std"] == 0
