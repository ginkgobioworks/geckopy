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

"""Integration layer with pytfa. Pytfa is not installed by default."""

import logging
import pickle
import zlib
from typing import Dict, Optional

import cobra
import pandas as pd
import pytfa
from pytfa.optim.variables import GenericVariable, LogConcentration

import geckopy


LOGGER = logging.getLogger(__name__)
PROT_DATA = {
    "pKa": [7],
    "deltaGf_err": 0,
    "mass_std": 0.0,
    "id": "protein",
    "nH_std": 12,
    "name": "protein",
    "formula": "H",
    "deltaGf_std": 0,
    "error": "Nil",
    "charge_std": 0,
    "struct_cues": {"ART_PROT": 0},
}
PROT_CUE_DATA = {
    "datfile": "ART.dgs",
    "error": 0,
    "formula": "",
    "charge": 0,
    "id": "ART_PROT",
    "small": True,
    "names": ["ART_PROT"],
    "energy": 0,
}


class ThermoModel(pytfa.ThermoModel):
    """Derived class to guard LC_vars against rewriting.

    This is required because we need to include `LC_vars` about proteins before
    calling `.convert()`, but convert reassigns it to an empty dictionary.
    """

    _LC_vars: Dict[cobra.Metabolite, GenericVariable] = {}

    @property
    def LC_vars(self):
        """Get dict LogConcentration variables to use in dGr constraints."""
        return self._LC_vars

    @LC_vars.setter
    def LC_vars(self, value: GenericVariable):
        if value:
            self._LC_vars = value


def adapt_gecko_to_thermo(
    ec_model: geckopy.Model,
    thermodb: Dict,
    compartment_data: Dict,
    solver: Optional[str] = None,
    *args,
    **kwargs,
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
    thermodb: Dict
        from pytfa.io.load_thermo. The format is explained at
        https://pytfa.readthedocs.io/en/latest/thermoDB.html
    compartment_data: Dict
        check https://pytfa.readthedocs.io/en/latest/model.html#compartment-data
    *args, **kwargs:
        which will be passed to the pytfa.ThermoModel.__init__
    """
    thermodb["metabolites"]["protein"] = PROT_DATA
    thermodb["cues"]["ART_PROT"] = PROT_CUE_DATA
    # preparing + converting the model sets the constrain linear coefficients
    # of proteins to usual values which may not be the intended behavior
    right_prot_coeffs = {
        prot.id: {
            var.name: coeff
            for var, coeff in ec_model.solver.constraints[prot.id]
            .get_linear_coefficients(ec_model.solver.constraints[prot.id].variables)
            .items()
        }
        for prot in ec_model.proteins
    }
    tmodel = ThermoModel(thermodb, ec_model, *args, **kwargs)
    if solver:
        tmodel.solver = solver
    tmodel.compartments = compartment_data
    # pseudometabolites of arm reactions are equaled to its opposite side
    for prot in tmodel.proteins:
        # pass ownership of the proteins to the ThermoModel
        prot._model = tmodel
        prot.thermo = pytfa.thermo.MetaboliteThermo(
            PROT_DATA,
            tmodel.compartments[prot.compartment]["pH"],
            tmodel.compartments[prot.compartment]["ionicStr"],
            tmodel.TEMPERATURE,
            tmodel.MIN_pH,
            tmodel.MAX_pH,
            tmodel.Debye_Huckel_B,
            tmodel.thermo_unit,
        )
    tmodel.prepare()
    # 0 dG formation for proteins so they are not used in thermo calculations
    # proteins are not included inside "model.metabolites" so we need to add'em
    for prot in tmodel.proteins:
        LC = tmodel.add_variable(LogConcentration, prot, lb=-1000, ub=1000)
        tmodel.LC_vars[prot] = LC
        prot.thermo.deltaGf_tr = 0
        prot.thermo.deltaGf_err = 0
    # restore whatever protein constraints we had previous to the model.repair
    for prot_id, constraint_terms in right_prot_coeffs.items():
        for var_name, coefficient in constraint_terms.items():
            tmodel.constraints[prot_id].set_linear_coefficients(
                {tmodel.variables[var_name]: coefficient}
            )
    tmodel.convert(verbose=False)
    # tmodel.convert(verbose=False, overwrite_lc_vars=False)
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
            met.annotation["seed_id"] = met.annotation["seed.id"]
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


def get_thermo_coverage(model: pytfa.thermo.ThermoModel, total=True):
    """Return the number of reactions that were assigned a thermodynamic variable.

    Parameters
    ----------
    model: pytfa.ThermoModel
    total: bool
        whether to report the total number of reactions covered or the percentage.
        Default: True
    """
    thermo_set = get_thermo_reactions(model)

    return len(thermo_set) if total else len(thermo_set) / len(model.reactions)


def get_thermo_reactions(model: pytfa.thermo.ThermoModel):
    r"""Return the set of reactions that were assigned a thermodynamic variable.

    In pytfa, a reaction will have thermodynamic constraints if and only if all
    of the metabolites in the reaction have the :math:`\Delta G_f`.
    """
    return {reac for reac in model.reactions if reac.thermo["computed"]}
