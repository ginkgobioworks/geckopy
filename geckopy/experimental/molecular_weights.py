"""Build dataframe with protein reactions identifiers, Uniprot IDs and MW."""
import re
import urllib.parse
import urllib.request
from typing import Dict, Optional

import pandas as pd

from geckopy.model import Model

__all__ = ["get_uniprot", "parse_mw", "extract_proteins"]


DEFAULT_PARAMS = {"from": "ACC+ID", "to": "ACC", "format": "txt", "query": ""}
URL = "https://www.uniprot.org/uploadlists/"
pat_mw = re.compile(r"\nSQ   SEQUENCE.+  (\d+) MW;")
pat_prot = re.compile(r"prot_(.+)")


def get_uniprot(query: str) -> str:
    """Get uniprot information corresponding to a query.

    Parameters
    ----------
    query: str
        an UNIPROT ID(s), separated by spaces

    """
    # WARNING: side effects on DEFAULT_PARAMS
    params = DEFAULT_PARAMS
    params["query"] = query
    data = urllib.parse.urlencode(params)
    data = data.encode("utf-8")
    req = urllib.request.Request(URL, data)
    with urllib.request.urlopen(req) as f:
        response = f.read()

    return response.decode("utf-8")


def parse_mw(uniprot_info: str) -> str:
    """Get all MW of uniprot text (Dalton)."""
    return pat_mw.findall(uniprot_info)


def _get_all_proteins(model: Model) -> Dict:
    """Generate a set of dict of Uniprot ID: reaction id from a `model`."""
    return {pat_prot.match(prot.id)[1]: prot.id for prot in model.proteins}


def extract_proteins(model, all_proteins: Optional[Dict] = None) -> pd.DataFrame:
    """Generate the dataframe protein reactions IDs, Uniprot IDs and MW.

    Parameters
    ----------
    model: cobra.Model
    all_proteins: dict
        dict of UNIPROT IDs to protein reaction identifiers as in the model.
        If None are supplied, the function will try to identify them with a
        simple regex.

    Returns
    -------
    df: pd.DataFrame

    """
    if all_proteins is None:
        all_proteins = _get_all_proteins(model)
    if not all_proteins:
        raise Exception("Set of proteins exchanges couldn't be resolved.")
    df = pd.DataFrame(
        {
            "uniprot": list(all_proteins.keys()),
            "reactions": list(all_proteins.values()),
        }
    )
    # get all the text in batch
    df["MW"] = parse_mw(get_uniprot(" ".join(list(all_proteins.keys()))))
    # from the regex, we get strings but we need numerics
    df["MW"] = pd.to_numeric(df["MW"])
    return df
