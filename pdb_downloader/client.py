"""RCSB search and file-download API."""

import re
from pathlib import Path

import requests

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
FILES_URL = "https://files.rcsb.org"
ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]*$")


def chemical_query(smiles: str, match: str = "sub-struct-graph-relaxed") -> dict:
    if not smiles.strip():
        raise ValueError("SMILES cannot be empty")
    if match not in {"graph-exact", "graph-relaxed", "sub-struct-graph-exact", "sub-struct-graph-relaxed"}:
        raise ValueError(f"Unsupported chemical match type: {match}")
    return {"type": "terminal", "service": "chemical", "parameters": {
        "value": smiles, "type": "descriptor", "descriptor_type": "SMILES", "match_type": match,
    }}


def metadata_query(protein: str | None = None, organism: str | None = None,
                   resolution: float | None = None) -> dict:
    nodes = []
    for attribute, value, operator in (
        ("rcsb_polymer_entity.pdbx_description", protein, "contains_phrase"),
        ("rcsb_entity_source_organism.scientific_name", organism, "exact_match"),
    ):
        if value is not None:
            if not value.strip():
                raise ValueError(f"{attribute} cannot be empty")
            nodes.append({"type": "terminal", "service": "text", "parameters": {
                "attribute": attribute, "operator": operator, "value": value,
            }})
    if resolution is not None:
        if not 0 < resolution < float("inf"):
            raise ValueError("Maximum resolution must be a positive finite number")
        nodes.append({"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_entry_info.resolution_combined",
            "operator": "less_or_equal", "value": resolution,
        }})
    if not nodes:
        raise ValueError("Provide at least one search criterion")
    return nodes[0] if len(nodes) == 1 else {"type": "group", "logical_operator": "and", "nodes": nodes}


class PDBClient:
    def __init__(self, session: requests.Session | None = None, timeout: float = 30):
        self.session = session or requests.Session()
        self.timeout = timeout

    def search(self, query: dict, return_type: str = "entry", *,
               page_size: int = 1000, max_results: int | None = None) -> list[str]:
        if return_type not in {"entry", "mol_definition"}:
            raise ValueError("return_type must be entry or mol_definition")
        if not 1 <= page_size <= 10000:
            raise ValueError("page_size must be between 1 and 10000")
        if max_results is not None and max_results < 1:
            raise ValueError("max_results must be positive")
        ids = []
        total = None
        while total is None or len(ids) < total:
            rows = min(page_size, max_results - len(ids)) if max_results is not None else page_size
            rows = min(rows, 10000 - len(ids))
            if rows <= 0:
                break
            payload = {"query": query, "return_type": return_type,
                       "request_options": {"paginate": {"start": len(ids), "rows": rows}}}
            response = self.session.post(SEARCH_URL, json=payload, timeout=self.timeout)
            if response.status_code == 204:
                if ids:
                    raise ValueError("RCSB returned an empty page before all results were retrieved")
                return []
            response.raise_for_status()
            data = response.json()
            batch = data.get("result_set", [])
            total = data["total_count"]
            if not batch and len(ids) < total:
                raise ValueError("RCSB returned an empty page before all results were retrieved")
            if not ids and total > 10000 and (max_results is None or max_results > 10000):
                payload["request_options"] = {"return_all_hits": True}
                response = self.session.post(SEARCH_URL, json=payload, timeout=self.timeout)
                response.raise_for_status()
                data = response.json()
                all_ids = [item["identifier"] for item in data.get("result_set", [])]
                if len(all_ids) != data["total_count"] or len(set(all_ids)) != len(all_ids):
                    raise ValueError("RCSB did not return a complete unique result set")
                return all_ids[:max_results] if max_results is not None else all_ids
            ids.extend(item["identifier"] for item in batch)
            if len(set(ids)) != len(ids):
                raise ValueError("RCSB returned duplicate identifiers across search pages")
        return ids

    def download(self, identifier: str, kind: str, file_format: str, output: Path,
                 *, overwrite: bool = False) -> Path:
        if not ID_PATTERN.fullmatch(identifier):
            raise ValueError(f"Invalid identifier: {identifier!r}")
        identifier = identifier.upper()
        if kind == "entry" and file_format in {"cif", "pdb"}:
            url = f"{FILES_URL}/download/{identifier}.{file_format}"
        elif kind == "component" and file_format in {"cif", "sdf"}:
            suffix = "_ideal.sdf" if file_format == "sdf" else ".cif"
            url = f"{FILES_URL}/ligands/download/{identifier}{suffix}"
        elif kind == "bird" and file_format == "cif" and identifier.startswith("PRD_"):
            url = f"{FILES_URL}/birds/download/{identifier}.cif"
        else:
            raise ValueError("Supported formats: entry cif/pdb; component cif/sdf; bird PRD_ identifier cif")
        output = Path(output)
        output.mkdir(parents=True, exist_ok=True)
        path = output / f"{identifier}.{file_format}"
        if path.exists() and not overwrite:
            raise FileExistsError(f"Output already exists: {path}")
        response = self.session.get(url, timeout=self.timeout)
        response.raise_for_status()
        if not response.content:
            raise ValueError(f"Empty response from {url}")
        path.write_bytes(response.content)
        return path
