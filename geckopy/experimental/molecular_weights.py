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

import gzip
import io
import re
from time import sleep
from typing import Callable, Dict, List, Optional

import pandas as pd
import requests
from tqdm import tqdm

from geckopy.model import Model


__all__ = ["get_uniprot", "parse_mw", "extract_proteins"]


pat_mw = re.compile(r"\nSQ   SEQUENCE.+  (\d+) MW;")
UNIPROT_PATTERN = re.compile(
    r"(?:prot_)?([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})"
)
protein_weights = {
    "A": 89.0932,
    "C": 121.1582,
    "D": 133.1027,
    "E": 147.1293,
    "F": 165.1891,
    "G": 75.0666,
    "H": 155.1546,
    "I": 131.1729,
    "K": 146.1876,
    "L": 131.1729,
    "M": 149.2113,
    "N": 132.1179,
    "O": 255.3134,
    "P": 115.1305,
    "Q": 146.1445,
    "R": 174.201,
    "S": 105.0926,
    "T": 119.1192,
    "U": 168.0532,
    "V": 117.1463,
    "W": 204.2252,
    "Y": 181.1885,
}


def get_uniprot(query: str) -> str:
    """Get uniprot information in JSON corresponding to a query."""
    url = f"https://rest.uniprot.org/uniprotkb/search?query={query}&format=json"
    return requests.get(url).json()["results"]


def parse_mw(uniprot_info: str) -> Dict[str, float]:
    """Get all MW of uniprot text (Dalton)."""
    return {
        entry["primaryAccession"]: float(entry["sequence"]["molWeight"])
        for entry in uniprot_info
    }


def _get_all_proteins(model: Model, key_fn) -> Dict[str, str]:
    """Generate a set of dict of Uniprot ID: reaction id from a `model`."""
    return {key_fn(prot): prot.id for prot in model.proteins}


def _molecular_weight(seq: str) -> float:
    """Calculate the molecular mass of DNA, RNA or protein sequences as float.

    Brought from biopython to avoid including the whole dependency. Only
    unambiguous letters are allowed.
    """
    seq = "".join(str(seq).split()).upper()  # Do the minimum formatting
    weight_table = protein_weights
    water = 18.010565

    try:
        weight = sum(weight_table[x] for x in seq) - (len(seq) - 1) * water
    except KeyError as e:
        raise ValueError(
            f"'{e}' is not a valid unambiguous letter for proteins"
        ) from None

    return weight


def download_uniprot_beta_data(
    accessions: List[str], fields: List[str]
) -> pd.DataFrame:
    # IF TOO MANY ACCESSIONS, RUN IN BATCHES
    # That's because running batches too big triggers errors.
    BATCH_SIZE = 2000
    if len(accessions) > BATCH_SIZE:
        accessions_batches = [
            accessions[i : i + BATCH_SIZE]
            for i in range(0, len(accessions), BATCH_SIZE)
        ]
        subdfs = [
            download_uniprot_beta_data(batch, fields)
            for batch in tqdm(accessions_batches)
        ]
        return pd.concat(subdfs, ignore_index=True)

    # SUBMIT THE QUERY

    data = {"ids": accessions, "from": "UniProtKB_AC-ID", "to": "UniProtKB"}
    job_submission_query = requests.post(
        "https://rest.uniprot.org/idmapping/run", data=data
    )
    job_id = job_submission_query.json()["jobId"]

    # POLL FOR STATUS
    endpoint = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        job_status_query = requests.get(endpoint)
        job_status = job_status_query.json()
        if "results" in job_status:
            break
        sleep(1)

    # GET THE RESULTS, PARSE THEM INTO A DATAFRAME
    endpoint = f"https://rest.uniprot.org/idmapping/uniprotkb/results/stream/{job_id}"
    response = requests.get(
        endpoint,
        params={
            "compressed": "true",
            "download": "true",
            "fields": ",".join(fields),
            "format": "tsv",
        },
    )
    gzip_reader = gzip.GzipFile(fileobj=io.BytesIO(response.content))
    return pd.read_csv(gzip_reader, sep="\t")


def extract_proteins(
    model,
    all_proteins: Optional[Dict] = None,
    key_fn: Callable[[str], str] = lambda x: UNIPROT_PATTERN.match(x.id)[1],
) -> pd.DataFrame:
    """Generate the dataframe protein reactions IDs, Uniprot IDs and MW.

    Parameters
    ----------
    model: cobra.Model)
    all_proteins: dict
        dict of UNIPROT IDs to protein reaction identifiers as in the model.
        If None are supplied, the function will try to identify them with a simple regex.
    key_fn: function
        mapping to extract the uniprot id from the protein. Default: regex
        matching protein id.

    Returns
    -------
    df: pd.DataFrame
        of columns `[uniprot, protein_id, MW, Sequence]`; where `protein_id` is
        the id in the model



    """
    if all_proteins is None:
        all_proteins = _get_all_proteins(model, key_fn)
    if not all_proteins:
        raise Exception("Set of proteins exchanges couldn't be resolved.")
    dfseqs = download_uniprot_beta_data(
        accessions=list(all_proteins.keys()), fields=["sequence"]
    )
    df_prot: pd.DataFrame = (
        pd.DataFrame(
            {"uniprot": all_proteins.keys(), "protein_id": all_proteins.values()}
        )
        .merge(dfseqs, how="outer", left_on="uniprot", right_on="From")
        .drop(columns=["From"])
    )

    df_prot["MW"] = 0
    for index, row in df_prot.iterrows():
        df_prot.loc[index, "MW"] = _molecular_weight(row["Sequence"])
    return df_prot[["uniprot", "protein_id", "MW", "Sequence"]]
