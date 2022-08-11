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

"""Build dataframe with protein reactions identifiers, Uniprot IDs and MW."""

import re
from itertools import zip_longest
from typing import Any, Dict, Iterable, Optional

import pandas as pd
import requests

from geckopy.model import Model


__all__ = ["get_uniprot", "parse_mw", "extract_proteins"]


pat_mw = re.compile(r"\nSQ   SEQUENCE.+  (\d+) MW;")
UNIPROT_PATTERN = re.compile(
    r"(?:prot_)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
)


def get_uniprot(query: str) -> str:
    """Get uniprot information in JSON corresponding to a query."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json"
    return requests.get(url).json()["results"]


def parse_mw(uniprot_info: str) -> dict[str, float]:
    """Get all MW of uniprot text (Dalton)."""
    return {
        entry["primaryAccession"]: float(entry["sequence"]["molWeight"])
        for entry in uniprot_info
    }


def _get_all_proteins(model: Model, key_fn) -> Dict:
    """Generate a set of dict of Uniprot ID: reaction id from a `model`."""
    return {key_fn(prot): prot.id for prot in model.proteins}


def _chunk_list(x: list[Any], n=25) -> list[Iterable[Any]]:
    return list(
        [el for el in chunk if el != ""]
        for chunk in zip_longest(*[iter(x)] * n, fillvalue="")
    )


def extract_proteins(
    model,
    all_proteins: Optional[Dict] = None,
    key_fn=lambda x: UNIPROT_PATTERN.match(x.id)[1],
) -> pd.DataFrame:
    """Generate the dataframe protein reactions IDs, Uniprot IDs and MW.

    Parameters
    ----------
    model: cobra.Model
    all_proteins: dict
        dict of UNIPROT IDs to protein reaction identifiers as in the model.
        If None are supplied, the function will try to identify them with a
        simple regex.
    key_fn: function
        mapping to extract the uniprot id from the protein. Default: regex
        matching on the protein id.

    Returns
    -------
    df: pd.DataFrame

    """
    if all_proteins is None:
        all_proteins = _get_all_proteins(model, key_fn)
    if not all_proteins:
        raise Exception("Set of proteins exchanges couldn't be resolved.")
    # get all the text in batches of 25 (max allowed by Uniprot)
    prot_to_mw = {}
    for chunk in _chunk_list(list(all_proteins.keys()), 25):
        prot_to_mw.update(
            parse_mw(get_uniprot("+OR+".join([f"accession:{prot}" for prot in chunk])))
        )
    df = pd.DataFrame(
        {
            "uniprot": list(all_proteins.keys()),
            "reactions": list(all_proteins.values()),
            "MW": [prot_to_mw[prot] for prot in all_proteins],
        }
    )
    # from the regex, we get strings but we need numerics
    df["MW"] = pd.to_numeric(df["MW"])
    return df
