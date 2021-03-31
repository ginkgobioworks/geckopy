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

"""Ensure that geckopy is consistent with cobrapy."""

import pandas as pd
import pytest

import geckopy


def test_unconstrained_ec_model_is_cobra_model(slim_solution, cobra_model):
    """Check that unconstrained ec_model returns the same maximum as the plain model."""
    assert round(cobra_model.slim_optimize(), 4) == round(slim_solution, 4)


@pytest.mark.xfail(reason="Loading cobrapy is not implemented")
def test_constrained_ec_model_is_not_cobra_model(cobra_model, experimental_copy_number):
    """Check that constrained ec_model returns different maximum than the plain model."""
    raw_proteomics = pd.read_csv(experimental_copy_number)
    ec_model = geckopy.experimental.from_copy_number(
        cobra_model.copy(),
        index=raw_proteomics["uniprot"],
        cell_copies=raw_proteomics["copies_per_cell"],
        stdev=raw_proteomics["stdev"],
        vol=2.3,
        dens=1.105e-12,
        water=0.3,
    )
    assert round(cobra_model.slim_optimize(), 4) != round(ec_model.slim_optimize(), 4)


@pytest.mark.xfail(reason="Loading cobrapy is not implemented")
def test_from_cobrapy_works(cobra_model):
    """Generate ec_gem from cobrapy_model."""
    ec_model = geckopy.Model(cobra_model)
    assert len(ec_model.proteins) == 1259
