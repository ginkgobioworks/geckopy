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
from geckopy.experimental.molecular_weights import extract_proteins
from geckopy.experimental.relaxation import (
    elastic_upper_relaxation,
    get_upper_relaxation,
    relax_proteomics_greedy,
)


def test_ec_model_from_copy_number_reduces_grow(
    ec_model_core, experimental_copy_number
):
    """Test that experimentally constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core,
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    sol = ec_model.slim_optimize()
    assert sol < 0.1 if not isnan(sol) else True


def test_model_from_mmol_gDW_cannot_grow(ec_model_core, experimental_mmol_gDW):
    """Test that experimentally constraining the model reduces the objective value."""
    processed_proteomics: pd.DataFrame = pd.read_csv(experimental_mmol_gDW)
    # proteins ID in the model are "prot_UNIPROT"
    processed_proteomics.index = processed_proteomics["uniprot"].apply(
        lambda x: f"prot_{x}"
    )
    processed_proteomics: pd.Series = processed_proteomics["mmol_per_cell"]
    ec_model = from_mmol_gDW(ec_model_core, processed_proteomics)
    objective_value = ec_model.slim_optimize()
    assert (objective_value < 0.1) if not isnan(objective_value) else True


def test_relaxed_ec_model_from_copy_number_can_grow(
    ec_model_core, experimental_copy_number
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    # proteins ID in the model are "prot_UNIPROT"
    ec_model = from_copy_number(
        ec_model_core,
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
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
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
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
    for prot_id in iis:
        ec_model.proteins.get_by_id(prot_id).concentration = None
    relaxed_sol = ec_model.slim_optimize()
    assert relaxed_sol >= 0.2 and relaxed_sol > sol


def test_extract_proteins_retrieves_all_mw(ec_model_core):
    """Test that all MW can be retrieved fom uniprot."""
    all_proteins = {prot.id[5:]: prot.id for prot in ec_model_core.proteins}
    df = extract_proteins(ec_model_core, all_proteins=all_proteins)
    assert len(df["MW"]) == len(ec_model_core.proteins)
    assert (df["MW"] != 0).all()
    for row in df.itertuples():
        ec_model_core.proteins.get_by_id(row[2]).mw = row[3]
