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
import pytest
import pytfa
from pytfa.io import load_thermoDB

from geckopy.experimental import from_copy_number
from geckopy.experimental.relaxation import Objective_rule
from geckopy.integration.relaxation import (
    relax_thermo_concentrations_proteins,
    relax_thermo_proteins,
)
from geckopy.integration.pytfa import adapt_gecko_to_thermo


def test_relax_thermo_dgr_and_proteins_works(
    ec_model_core,
    experimental_copy_number,
    thermodb,
    compartment_data,
    slim_solution_core,
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core.copy(),
        # proteins in the model are in prot_UNIPROT form
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    thermodb = load_thermoDB(thermodb)
    compartment_data = pytfa.io.read_compartment_data(compartment_data)
    tmodel = adapt_gecko_to_thermo(ec_model, thermodb, compartment_data)
    sol_pre_relax = tmodel.slim_optimize()
    tmodel.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = slim_solution_core
    # this is an inplace operation so we have to regenerate the model after it
    iis_sum_obj, status = relax_thermo_proteins(
        tmodel,
        [
            prot.id
            for prot in tmodel.proteins
            if prot.concentration is not None and not isnan(prot.concentration)
        ],
        Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
    )
    assert status == "optimal"
    ec_model = from_copy_number(
        ec_model_core.copy(),
        # proteins in the model are in prot_UNIPROT form
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    tmodel = adapt_gecko_to_thermo(ec_model, thermodb, compartment_data)
    for var in iis_sum_obj:
        if var in tmodel.proteins:
            tmodel.proteins.get_by_id(var).concentration = None
    sol_post_relax = tmodel.slim_optimize()
    assert sol_pre_relax < 0.1 if not isnan(sol_pre_relax) else True
    assert sol_post_relax > 0.1


def test_relax_concentrations_and_proteins_works(
    ec_model_core,
    experimental_copy_number,
    thermodb,
    compartment_data,
    slim_solution_core,
):
    """Test that constraining the model works."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = from_copy_number(
        ec_model_core.copy(),
        # proteins in the model are in prot_UNIPROT form
        index=raw_proteomics["uniprot"].apply(lambda x: f"prot_{x}"),
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    thermodb = load_thermoDB(thermodb)
    compartment_data = pytfa.io.read_compartment_data(compartment_data)
    tmodel = adapt_gecko_to_thermo(ec_model, thermodb, compartment_data)
    sol_pre_relax = tmodel.slim_optimize()
    tmodel.variables.LC_atp_c.ub = 2e2
    tmodel.variables.LC_atp_c.lb = 1e2
    tmodel.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = slim_solution_core * 0.4
    # this is an inplace operation so we have to regenerate the model after it
    iis_sum_obj, status = relax_thermo_concentrations_proteins(
        tmodel,
        [
            prot.id
            for prot in tmodel.proteins
            if prot.concentration is not None and not isnan(prot.concentration)
        ],
        Objective_rule.MIN_ELASTIC_SUM_OBJECTIVE,
    )
    tmodel.objective = tmodel.reactions.BIOMASS_Ecoli_core_w_GAM
    sol_post_relax = tmodel.slim_optimize()
    # in the MILP relaxation atp is relaxed, whereas in the SUM LP OBJ relaxation
    # adp, gln and pi gets relaxed (in the opposite direction) compensating atp
    assert "NegSlackLC_atp_c" in iis_sum_obj or "PosSlackLC_adp_c" in iis_sum_obj
    assert status == "optimal"
    assert sol_pre_relax < 0.1 if not isnan(sol_pre_relax) else True
    assert pytest.approx(sol_post_relax) == slim_solution_core * 0.4
