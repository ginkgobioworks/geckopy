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

"""Loading experimental data."""

import logging
from typing import Union

import cobra
import pandas as pd

from geckopy.model import Model


__all__ = ["limit_proteins", "from_mmol_gDW", "from_copy_number"]
LOGGER = logging.getLogger(__name__)


def limit_proteins(model: Union[Model, cobra.Model], measurements: pd.Series):
    """Apply proteomics `measurements` to `model` (in-place).

    Adapted from
    https://github.com/DD-DeCaF/simulations/blob/devel/src/simulations/modeling/driven.py

    Parameters
    ----------
    model: cobra.Model (or geckopy.Model)
        The enzyme-constrained model.
    measurements: pd.Series
        Protein abundances in mmol / gDW.
        - Index: protein ID as it appears in the models.
        - Values: protein abundance in mmol / gDW

    Returns
    -------
    applied_measurements: int
        number of proteins that were effectively constrained in the model

    Warnings
    --------
    A warning will be logged If a protein in the measurements cannot be found
    in the model.

    """
    is_ec_model = hasattr(model, "proteins")
    applied_measurements = 0
    for protein_id, measure in measurements.items():
        try:
            rxn = (
                model.proteins.get_by_id(protein_id)
                if is_ec_model
                else model.reactions.get_by_id(f"{protein_id}_exchange")
            )
        except KeyError:
            LOGGER.warn(
                f"protein {protein_id} from measurements was not found in the model"
            )
        else:
            # update only upper_bound (as enzymes can be unsaturated):
            if is_ec_model:
                rxn.concentration = measure
            else:
                rxn.upper_bound = measure
            applied_measurements += 1
    return applied_measurements


def from_mmol_gDW(model: cobra.Model, processed_proteomics: pd.Series) -> cobra.Model:
    """Apply proteomics constraints to model and return the copied EC model.

    Parameters
    ----------
    measurements: pd.Series
        Protein abundances in mmol / gDW.
        - Index: protein ID as it appears in the models.
        - Values: protein abundance in mmol / gDW
    """
    ec_model = Model(model.copy())
    _ = limit_proteins(ec_model, processed_proteomics)
    return ec_model


def from_copy_number(
    model: cobra.Model,
    index: pd.Series,
    cell_copies: pd.Series,
    stdev: pd.Series,
    vol: float,
    dens: float,
    water: float,
) -> cobra.Model:
    """Convert `cell_copies` to mmol/gDW and apply them to `model`.

    Parameters
    ----------
    model: cobra.Model
        cobra or geckopy Model (will be converted to geckopy.Model). It is NOT
        modified inplace.
    index: pd.Series
        uniprot IDs
    cell_copies: pd.Series
        cell copies/ cell per proteins
    stdev: pd.Series
        standard deviation of the cell copies
    vol: float
        cell volume
    dens: float
        cell density
    water: float
        water content fraction (0-1)

    Returns
    -------
    geckopy.Model
        with the proteomics constraints applied
    """
    df = pd.DataFrame({"cell_copies": cell_copies, "CV": stdev})
    # from molecules/cell to mmol/gDW
    df["copies_upper"] = df["cell_copies"] + 0.5 * df["CV"] / 100 * df["cell_copies"]
    df["mmol_per_cell"] = df["copies_upper"] / 6.022e21
    proteomics = df["mmol_per_cell"] / (vol * dens * water)
    proteomics.index = index
    return from_mmol_gDW(model, proteomics)
