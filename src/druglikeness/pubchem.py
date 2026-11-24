"""Retrieve ligand structures from PubChem PUG REST."""

import csv
import json
import time
import urllib.error
import urllib.request
from dataclasses import dataclass
from pathlib import Path

PROPERTY_FIELDS = (
    "Title",
    "IUPACName",
    "SMILES",
    "MolecularFormula",
    "MolecularWeight",
    "XLogP",
    "TPSA",
    "HBondDonorCount",
    "HBondAcceptorCount",
    "RotatableBondCount",
)

PUG_REST = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid"
USER_AGENT = (
    "ligand-druglikeness/0.2 "
    "(https://github.com/Ishitta-Sarkar/computational-drug-discovery-platform)"
)


@dataclass(frozen=True)
class LigandQuery:
    name: str
    pubchem_cid: int


@dataclass(frozen=True)
class ReferenceProperties:
    """PubChem's own property values. They are not inputs to the filters."""

    molecular_weight: float
    xlogp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int


@dataclass(frozen=True)
class CompoundRecord:
    name: str
    pubchem_cid: int
    pubchem_title: str
    iupac_name: str
    smiles: str
    molecular_formula: str
    reference: ReferenceProperties


@dataclass(frozen=True)
class Retrieval:
    retrieved_utc: str
    request_url: str
    raw: dict
    records: tuple[CompoundRecord, ...]


def load_queries(path: Path) -> list[LigandQuery]:
    queries: list[LigandQuery] = []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        if "docking_score" in fieldnames:
            raise ValueError(
                "data/ligands.csv contains docking_score. A docking score belongs "
                "to one ligand-target calculation. This workflow does not dock, "
                "and it will not read an unattached score."
            )
        required = {"compound_name", "pubchem_cid"}
        missing = required.difference(fieldnames)
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for line_number, row in enumerate(reader, start=2):
            name = (row.get("compound_name") or "").strip()
            cid_text = (row.get("pubchem_cid") or "").strip()
            if not name or not cid_text:
                raise ValueError(f"{path} line {line_number} is missing a name or CID.")
            if not cid_text.isdigit():
                raise ValueError(f"{path} line {line_number} has a non-integer CID: {cid_text}")
            queries.append(LigandQuery(name=name, pubchem_cid=int(cid_text)))
    if not queries:
        raise ValueError(f"{path} has no ligands.")
    cids = [query.pubchem_cid for query in queries]
    if len(cids) != len(set(cids)):
        raise ValueError(f"{path} repeats a PubChem CID.")
    return queries


def build_url(cids: list[int] | tuple[int, ...]) -> str:
    cid_path = ",".join(str(cid) for cid in cids)
    properties = ",".join(PROPERTY_FIELDS)
    return f"{PUG_REST}/{cid_path}/property/{properties}/JSON"


def fetch_records(queries: list[LigandQuery]) -> Retrieval:
    from datetime import datetime, timezone

    url = build_url([query.pubchem_cid for query in queries])
    raw = _get_json(url)
    records = parse_property_table(raw, queries)
    retrieved = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    return Retrieval(
        retrieved_utc=retrieved,
        request_url=url,
        raw=raw,
        records=tuple(records),
    )


def parse_property_table(raw: dict, queries: list[LigandQuery]) -> list[CompoundRecord]:
    try:
        rows = raw["PropertyTable"]["Properties"]
    except (KeyError, TypeError) as exc:
        raise ValueError("PubChem JSON has no PropertyTable.Properties list.") from exc
    if not isinstance(rows, list):
        raise ValueError("PubChem PropertyTable.Properties is not a list.")

    by_cid: dict[int, dict] = {}
    for row in rows:
        cid = int(row["CID"])
        if cid in by_cid:
            raise ValueError(f"PubChem returned CID {cid} more than once.")
        by_cid[cid] = row

    records: list[CompoundRecord] = []
    for query in queries:
        row = by_cid.get(query.pubchem_cid)
        if row is None:
            raise ValueError(f"PubChem response has no CID {query.pubchem_cid} ({query.name}).")
        title = str(row.get("Title", "")).strip()
        if title.casefold() != query.name.casefold():
            raise ValueError(
                f"CID {query.pubchem_cid} is titled {title!r} at PubChem, "
                f"not {query.name!r}. Refusing to attach this structure."
            )
        smiles = str(row.get("SMILES", "")).strip()
        formula = str(row.get("MolecularFormula", "")).strip()
        iupac = str(row.get("IUPACName", "")).strip()
        if not smiles or not formula:
            raise ValueError(f"CID {query.pubchem_cid} is missing SMILES or a formula.")
        reference = ReferenceProperties(
            molecular_weight=float(row["MolecularWeight"]),
            xlogp=float(row["XLogP"]),
            tpsa=float(row["TPSA"]),
            hbd=int(row["HBondDonorCount"]),
            hba=int(row["HBondAcceptorCount"]),
            rotatable_bonds=int(row["RotatableBondCount"]),
        )
        records.append(
            CompoundRecord(
                name=query.name,
                pubchem_cid=query.pubchem_cid,
                pubchem_title=title,
                iupac_name=iupac,
                smiles=smiles,
                molecular_formula=formula,
                reference=reference,
            )
        )
    return records


def _get_json(url: str, attempts: int = 3) -> dict:
    delay_seconds = 2.0
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(
            url,
            headers={"Accept": "application/json", "User-Agent": USER_AGENT},
        )
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                payload = response.read().decode("utf-8")
            return json.loads(payload)
        except urllib.error.HTTPError as exc:
            last_error = exc
            body = exc.read().decode("utf-8", errors="replace")[:500]
            retriable = exc.code in {429, 500, 502, 503, 504}
            if not retriable or attempt == attempts - 1:
                raise RuntimeError(
                    f"PubChem request failed with HTTP {exc.code} for {url}: {body}"
                ) from exc
        except urllib.error.URLError as exc:
            last_error = exc
            if attempt == attempts - 1:
                raise RuntimeError(f"PubChem request failed for {url}: {exc.reason}") from exc
        time.sleep(delay_seconds)
        delay_seconds *= 2
    raise RuntimeError(f"PubChem request failed for {url}") from last_error
