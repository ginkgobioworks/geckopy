"""Integration layer with pytfa. Pytfa is not installed by default."""
import pickle
import zlib
from typing import Dict

import cobra
import pandas as pd

import geckopy


def prepare_gecko_to_thermo(
    model: geckopy.Model,
    thermodb: Dict,
    mnx_file: str,
):
    """Prepare the annotation of metabolites and proteins of `model` to pytfa.

    It modifies the `model` and the `thermodb` in place.

    The 'seed_id' annotation is hardcoded in pytfa so the metabolites and the
    database must have that annotation to be used. 'metanetx.chemical' annotation
    are translated to seed if possible.

    If proteins were added without dG energy, they assumed to be missing metabolites
    so the reactions with enzymes would be ignored. Adding 0 formation energy
    results in ignoring proteins for dG, since catalyzers do not affect dG.

    Parameters
    ----------
    model: geckopy.Model
    thermodb: dict
        from pytfa.io.load_thermo. The format is explained at
        https://pytfa.readthedocs.io/en/latest/thermoDB.html
    mnx_file: str
        path to "chem_xref.tsv" file. It can be downloaded from
        https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv
    """
    translate_model_mnx_to_seed(model, thermodb, mnx_file)
    add_dummy_protein_info(model, thermodb)


def translate_model_mnx_to_seed(
    model: cobra.Model,
    thermodb: Dict,
    mnx_file: str,
):
    """Add a seed_id annotation to every metabolite."""
    df = pd.read_csv(
        mnx_file,
        sep="\t",
        comment="#",
        names=["xref", "mnx", "annotation"],
    )
    # perf filter
    df = df.loc[df.xref.str.startswith("seed"), :]
    df.xref = df.xref.apply(lambda x: x[-8:])
    for met in model.metabolites:
        if "seed.id" in met.annotation:
            met.annotation["seed.id"] = met.annotation["seed_id"]
        if "seed_id" in met.annotation:
            continue
        if "metanetx.chemical" in met.annotation:
            mnx_id = met.annotation["metanetx.chemical"]
            if isinstance(mnx_id, str):
                annotations = df.loc[
                    df.mnx == met.annotation["metanetx.chemical"], "xref"
                ]
            else:
                annotations = df.loc[
                    df.mnx.isin(met.annotation["metanetx.chemical"]), "xref"
                ]
            annotation = [ann for ann in annotations if ann in thermodb["metabolites"]]
            if annotation:
                # just pick the first one
                met.annotation["seed_id"] = annotation[0]


def add_dummy_protein_info(model: geckopy.Model, thermodb: Dict):
    """Add 0 formation energy for proteins."""
    for prot in model.proteins:
        prot.annotation["seed_id"] = "protein"
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
    }


def write_thermodb(thermodb: Dict, filename: str):
    """Deserialize the `thermoDB` to a compressed file at `filename`."""
    with open(filename, "wb") as f:
        f.write(zlib.compress(pickle.dumps(thermodb)))
