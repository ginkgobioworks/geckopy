"""Integration layer with pytfa. Pytfa is not installed by default."""
import pickle
import zlib
from typing import Dict, Optional

import cobra
import pandas as pd
import pytfa

import geckopy


def adapt_gecko_to_thermo(
    ec_model: geckopy.Model,
    thermodb: dict,
    compartment_data: dict,
    solver: Optional[str] = None,
    *args,
    **kwargs
) -> pytfa.ThermoModel:
    """Prepare and convert gecko model to `pytfa.ThermoModel`.

    It modifies the `model` and the `thermodb` in place.

    The 'seed_id' annotation is hardcoded in pytfa so the metabolites and the
    database must have that annotation to be used.

    If proteins were added without dG energy, they would be treated as missing
    metabolites so the reactions with enzymes would be ignored. Adding 0
    formation energy results in ignoring proteins for dG, since catalyzers do
    not affect dG.

    Parameters
    ----------
    model: geckopy.Model
    thermodb: dict
        from pytfa.io.load_thermo. The format is explained at
        https://pytfa.readthedocs.io/en/latest/thermoDB.html
    compartment_data: dict
        check https://pytfa.readthedocs.io/en/latest/model.html#compartment-data
    *args, **kwargs:
        which will be passed to the pytfa.ThermoModel.__init__
    """
    tmodel = pytfa.ThermoModel(thermodb, ec_model, *args, **kwargs)
    if solver:
        tmodel.solver = solver
    tmodel.compartments = compartment_data
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
    prot_data = thermodb["metabolites"]["protein"]
    for prot in tmodel.proteins:
        CompartmentpH = tmodel.compartments[prot.compartment]["pH"]
        CompartmentionicStr = tmodel.compartments[prot.compartment]["ionicStr"]
        prot.thermo = pytfa.thermo.MetaboliteThermo(
            prot_data,
            CompartmentpH,
            CompartmentionicStr,
            tmodel.TEMPERATURE,
            tmodel.MIN_pH,
            tmodel.MAX_pH,
            tmodel.Debye_Huckel_B,
            tmodel.thermo_unit,
        )
    tmodel.prepare()
    for prot in tmodel.proteins:
        # metabolites with this formula are ignored
        prot.formula = "H"
    tmodel.convert(verbose=False)
    return tmodel


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


def write_thermodb(thermodb: Dict, filename: str):
    """Deserialize the `thermoDB` to a compressed file at `filename`."""
    with open(filename, "wb") as f:
        f.write(zlib.compress(pickle.dumps(thermodb)))
