"""Read the ligand list and fetch or load PubChem compound records.

The ligand file stores a name, a PubChem CID and a chemical class. SMILES
and the other structure fields come from PubChem PUG REST, either live or
from the JSON cache written by ``scripts/fetch_pubchem.py``.
"""

from __future__ import annotations

import csv
import json
import urllib.error
import urllib.request
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

PUG_PROPERTIES = (
    "Title",
    "IUPACName",
    "MolecularFormula",
    "InChI",
    "InChIKey",
    "SMILES",
    "ConnectivitySMILES",
)

USER_AGENT = (
    "ligand-druglikeness-screen/0.2 (educational cheminformatics; mailto:ishittas4@gmail.com)"
)

SCHEMA_VERSION = 1


@dataclass(frozen=True)
class LigandQuery:
    """One row of the local ligand list, before any descriptor is computed."""

    compound_name: str
    pubchem_cid: int
    compound_class: str


@dataclass(frozen=True)
class PubChemCompound:
    """A PubChem record aligned to one ligand-list row."""

    compound_name: str
    compound_class: str
    pubchem_cid: int
    title: str
    iupac_name: str
    molecular_formula: str
    isomeric_smiles: str
    connectivity_smiles: str
    inchi: str
    inchikey: str


def load_library(path: Path) -> list[LigandQuery]:
    """Load ``compound_name,pubchem_cid,compound_class`` rows."""
    queries: list[LigandQuery] = []
    seen: set[int] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"compound_name", "pubchem_cid", "compound_class"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError(
                f"{path} must have columns compound_name, pubchem_cid, compound_class."
            )
        if "docking_score" in reader.fieldnames:
            raise ValueError(
                f"{path} contains a docking_score column. "
                "This screen does not accept docking scores."
            )
        for line_number, row in enumerate(reader, start=2):
            name = (row.get("compound_name") or "").strip()
            cid_text = (row.get("pubchem_cid") or "").strip()
            compound_class = (row.get("compound_class") or "").strip()
            if not name or not cid_text or not compound_class:
                raise ValueError(f"{path}:{line_number} has an empty required field.")
            if not cid_text.isdigit():
                raise ValueError(f"{path}:{line_number} has a non-integer CID: {cid_text}")
            cid = int(cid_text)
            if cid in seen:
                raise ValueError(f"{path} lists CID {cid} more than once.")
            seen.add(cid)
            queries.append(
                LigandQuery(
                    compound_name=name,
                    pubchem_cid=cid,
                    compound_class=compound_class,
                )
            )
    if not queries:
        raise ValueError(f"{path} does not list any ligands.")
    return queries


def pug_property_url(cids: list[int]) -> str:
    """PUG REST URL for the properties this screen needs."""
    cid_path = ",".join(str(cid) for cid in cids)
    properties = ",".join(PUG_PROPERTIES)
    return (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid_path}/property/{properties}/JSON"
    )


def fetch_property_table(
    cids: list[int],
    opener: Any = None,
    timeout: float = 60,
) -> tuple[list[dict[str, Any]], str]:
    """Request PubChem properties. ``opener`` replaces ``urlopen`` in tests."""
    if not cids:
        raise ValueError("No PubChem CIDs to fetch.")
    url = pug_property_url(cids)
    request = urllib.request.Request(
        url,
        headers={"User-Agent": USER_AGENT, "Accept": "application/json"},
    )
    open_fn = opener or urllib.request.urlopen
    try:
        with open_fn(request, timeout=timeout) as response:
            payload = json.load(response)
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")[:500]
        raise RuntimeError(f"PubChem request failed ({exc.code}) for {url}: {detail}") from exc
    except urllib.error.URLError as exc:
        raise RuntimeError(f"PubChem request failed for {url}: {exc.reason}") from exc
    try:
        properties = payload["PropertyTable"]["Properties"]
    except (KeyError, TypeError) as exc:
        raise RuntimeError(f"PubChem response did not contain a property table: {url}") from exc
    if not isinstance(properties, list):
        raise RuntimeError(f"PubChem property table was not a list: {url}")
    return properties, url


def _require_text(row: dict[str, Any], key: str, cid: int) -> str:
    value = row.get(key)
    if value is None or str(value).strip() == "":
        raise ValueError(f"PubChem CID {cid} is missing {key}.")
    return str(value).strip()


def records_from_property_table(
    queries: list[LigandQuery],
    properties: list[dict[str, Any]],
) -> list[PubChemCompound]:
    """Align PubChem rows to the ligand list and check each title."""
    by_cid = {int(row["CID"]): row for row in properties if "CID" in row}
    compounds: list[PubChemCompound] = []
    for query in queries:
        row = by_cid.get(query.pubchem_cid)
        if row is None:
            raise ValueError(
                f"PubChem returned no record for CID {query.pubchem_cid} ({query.compound_name})."
            )
        title = _require_text(row, "Title", query.pubchem_cid)
        if title.casefold() != query.compound_name.casefold():
            raise ValueError(
                f"PubChem title for CID {query.pubchem_cid} is {title!r}, "
                f"which does not match the library name {query.compound_name!r}."
            )
        compounds.append(
            PubChemCompound(
                compound_name=query.compound_name,
                compound_class=query.compound_class,
                pubchem_cid=query.pubchem_cid,
                title=title,
                iupac_name=_require_text(row, "IUPACName", query.pubchem_cid),
                molecular_formula=_require_text(row, "MolecularFormula", query.pubchem_cid),
                isomeric_smiles=_require_text(row, "SMILES", query.pubchem_cid),
                connectivity_smiles=_require_text(row, "ConnectivitySMILES", query.pubchem_cid),
                inchi=_require_text(row, "InChI", query.pubchem_cid),
                inchikey=_require_text(row, "InChIKey", query.pubchem_cid),
            )
        )
    return compounds


def fetch_compounds(
    queries: list[LigandQuery],
    opener: Any = None,
    timeout: float = 60,
) -> tuple[list[PubChemCompound], str]:
    """Fetch every CID in ``queries`` and return compounds in library order."""
    properties, url = fetch_property_table(
        [query.pubchem_cid for query in queries],
        opener=opener,
        timeout=timeout,
    )
    return records_from_property_table(queries, properties), url


def cache_document(
    compounds: list[PubChemCompound],
    request_url: str,
    retrieved_on: str,
) -> dict[str, Any]:
    """JSON-serialisable cache with the request that produced it."""
    return {
        "schema_version": SCHEMA_VERSION,
        "source": "PubChem PUG REST",
        "retrieved_on": retrieved_on,
        "request_url": request_url,
        "compounds": [asdict(compound) for compound in compounds],
    }


def write_cache(path: Path, document: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")


def load_cache(path: Path) -> tuple[dict[str, Any], list[PubChemCompound]]:
    """Load a cache written by ``write_cache``."""
    try:
        document = json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"PubChem cache not found at {path}. Run python scripts/fetch_pubchem.py."
        ) from exc
    if document.get("schema_version") != SCHEMA_VERSION:
        raise ValueError(f"{path} has an unsupported schema_version.")
    compounds = [
        PubChemCompound(
            compound_name=row["compound_name"],
            compound_class=row["compound_class"],
            pubchem_cid=int(row["pubchem_cid"]),
            title=row["title"],
            iupac_name=row["iupac_name"],
            molecular_formula=row["molecular_formula"],
            isomeric_smiles=row["isomeric_smiles"],
            connectivity_smiles=row["connectivity_smiles"],
            inchi=row["inchi"],
            inchikey=row["inchikey"],
        )
        for row in document["compounds"]
    ]
    return document, compounds


def assert_cache_matches_library(
    compounds: list[PubChemCompound],
    queries: list[LigandQuery],
) -> None:
    """Require the cache to list the same names and CIDs, in the same order."""
    cached = [(compound.compound_name, compound.pubchem_cid) for compound in compounds]
    wanted = [(query.compound_name, query.pubchem_cid) for query in queries]
    if cached != wanted:
        raise ValueError(
            "The PubChem cache does not match the ligand list. "
            "Run python scripts/fetch_pubchem.py to refresh it."
        )
