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

"""Tests for experimental data loading."""

from math import isnan

import pandas as pd

from geckopy.experimental import from_copy_number, from_mmol_gDW
from geckopy.experimental.relaxation import (
    elastic_upper_relaxation,
    get_upper_relaxation,
    relax_proteomics_greedy,
)


def test_ec_model_from_copy_number_reduces_grow(
    ec_model_core, experimental_copy_number
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    sol = ec_model.slim_optimize()
    assert sol < 0.1 if not isnan(sol) else True


def test_model_from_mmol_gDW_can_grow(ec_model_core, experimental_mmol_gDW):
    """Test that constraining the model works."""
    processed_proteomics = pd.read_csv(experimental_mmol_gDW)
    ec_model = from_mmol_gDW(ec_model_core, processed_proteomics)
    assert ec_model.slim_optimize() > 0.0


def test_relaxed_ec_model_from_copy_number_can_grow(
    ec_model_core, experimental_copy_number
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    sol = ec_model.slim_optimize()
    sol = sol if not isnan(sol) and sol else 0.0
    ec_model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.8
    iis, _ = get_upper_relaxation(
        ec_model.copy(), [prot.id for prot in ec_model.proteins]
    )
    ec_model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.0
    final_gr, prots = relax_proteomics_greedy(ec_model, 0.5, protein_set=iis)
    relaxed_sol = ec_model.slim_optimize()
    assert relaxed_sol >= 0.5 and relaxed_sol > sol


def test_irreductibly_relaxed_ec_model_from_copy_number_can_grow(
    ec_model_core, experimental_copy_number
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core,
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    sol = ec_model.slim_optimize()
    sol = sol if not isnan(sol) and sol else 0.0
    ec_model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.8
    iis = elastic_upper_relaxation(ec_model, [prot.id for prot in ec_model.proteins])
    ec_model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.0
    final_gr, prots = relax_proteomics_greedy(ec_model, 0.2, protein_set=iis)
    relaxed_sol = ec_model.slim_optimize()
    assert relaxed_sol >= 0.2 and relaxed_sol > sol
