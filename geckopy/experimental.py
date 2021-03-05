"""Loading experimental data."""

import cobra
import pandas as pd

from .model import Model


def limit_proteins(model: cobra.Model, measurements: pd.DataFrame):
    """Apply proteomics `measurements` to `model`.

    Adapted from
    https://github.com/DD-DeCaF/simulations/blob/devel/src/simulations/modeling/driven.py

    Parameters
    ----------
    model: cobra.Model (or geckopy.Model)
        The enzyme-constrained model.
    measurements : pd.DataFrame
        Protein abundances in mmol / gDW.

    """
    is_ec_model = hasattr(model, "proteins")
    for protein_id, measure in measurements.items():
        try:
            rxn = (
                model.proteins.get_by_id(f"prot_{protein_id}")
                if is_ec_model
                else model.reactions.get_by_id(f"prot_{protein_id}_exchange")
            )
        except KeyError:
            pass
        else:
            # update only upper_bound (as enzymes can be unsaturated):
            if is_ec_model:
                rxn.concentration = measure
            else:
                rxn.upper_bound = measure


def from_mmol_gDW(
    model: cobra.Model, processed_proteomics: pd.DataFrame
) -> cobra.Model:
    """Apply proteomics constraints to model and return the copied EC model."""
    ec_model = Model(model)
    limit_proteins(ec_model, processed_proteomics)
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
    """Convert `cell_copies` to mmol/gDW and apply them to `model`."""
    df = pd.DataFrame({"cell_copies": cell_copies, "CV": stdev})
    # from molecules/cell to mmol/gDW
    df["copies_upper"] = df["cell_copies"] + 0.5 * df["CV"] / 100 * df["cell_copies"]
    df["mmol_per_cell"] = df["copies_upper"] / 6.022e21
    proteomics = df["mmol_per_cell"] / (vol * dens * water)
    proteomics.index = index
    return from_mmol_gDW(model, proteomics)
