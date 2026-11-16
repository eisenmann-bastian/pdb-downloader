"""Extract non-polymer occurrences with experimental coordinates and CCD bonds."""

import csv
from collections import defaultdict
from pathlib import Path

import gemmi
import requests
from rdkit import Chem

from .client import PDBClient

FIELDS = ["entry_id", "model", "label_asym_id", "auth_asym_id", "auth_seq_id",
          "insertion_code", "component_id", "alt_id", "file", "status", "detail"]
ORDERS = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
          "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}


def _column(category: dict, name: str, index: int, default: str = "") -> str:
    values = category.get(name)
    if values is None:
        return default
    value = values[index]
    return str(value) if value is not None and value is not False else default


def _component(path: Path) -> tuple[dict, list]:
    block = gemmi.cif.read_file(str(path)).sole_block()
    atoms = block.get_mmcif_category("_chem_comp_atom.")
    bonds = block.get_mmcif_category("_chem_comp_bond.")
    if not atoms or "charge" not in atoms:
        raise ValueError(f"Component file lacks atom definitions or formal charges: {path}")
    definitions = {}
    for i, atom_id in enumerate(atoms["atom_id"]):
        charge = _column(atoms, "charge", i)
        definitions[atom_id] = (_column(atoms, "type_symbol", i), int(charge) if charge else None)
    bond_list = []
    for i, atom_id in enumerate(bonds.get("atom_id_1", [])):
        bond_list.append((atom_id, _column(bonds, "atom_id_2", i),
                          _column(bonds, "value_order", i).upper()))
    return definitions, bond_list


def _molecule(rows: list[dict], definitions: dict, bonds: list, *,
              allow_unknown_charge: bool) -> tuple[Chem.Mol, list[str], list[str]]:
    molecule = Chem.RWMol()
    indices = {}
    unknown_charges = []
    conformer = Chem.Conformer(len(rows))
    for i, row in enumerate(rows):
        atom_id = row["atom_id"]
        if atom_id in indices:
            raise ValueError(f"Duplicate atom {atom_id} in a conformer")
        if atom_id not in definitions:
            raise ValueError(f"Atom {atom_id} is absent from the component definition")
        element, charge = definitions[atom_id]
        if charge is None:
            if not allow_unknown_charge:
                raise ValueError(f"Unknown CCD formal charge for {atom_id}; use --allow-unknown-charge to assume zero")
            unknown_charges.append(atom_id)
            charge = 0
        if element.upper() != row["element"].upper():
            raise ValueError(f"Element mismatch for atom {atom_id}: {row['element']} vs {element}")
        atom = Chem.Atom(element.capitalize())
        atom.SetFormalCharge(charge)
        atom.SetProp("_TriposAtomName", atom_id)
        indices[atom_id] = molecule.AddAtom(atom)
        conformer.SetAtomPosition(i, (float(row["x"]), float(row["y"]), float(row["z"])))
    for first, second, order in bonds:
        if first in indices and second in indices:
            if order not in ORDERS:
                raise ValueError(f"Unsupported CCD bond order {order}")
            molecule.AddBond(indices[first], indices[second], ORDERS[order])
    if len(rows) > 1 and molecule.GetNumBonds() == 0:
        raise ValueError("No CCD bonds connect the observed atoms")
    molecule.AddConformer(conformer)
    result = molecule.GetMol()
    Chem.SanitizeMol(result)
    missing = sorted(name for name, (element, _) in definitions.items()
                     if element.upper() not in {"H", "D"} and name not in indices)
    return result, missing, unknown_charges


def extract_ligands(cif: Path, output: Path, *, components: set[str] | None = None,
                    smarts: str | None = None, ccd_dir: Path | None = None,
                    client: PDBClient | None = None, deduplicate: bool = False,
                    allow_unknown_charge: bool = False) -> list[dict]:
    """Extract non-polymer residues; report every selected alternate occurrence.

    CCD CIF files are read from ccd_dir, or fetched into output/ccd when omitted.
    Failed occurrences are recorded in manifest.csv and cause a RuntimeError.
    """
    pattern = Chem.MolFromSmarts(smarts) if smarts is not None else None
    if smarts is not None and pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError(f"Extraction output directory must be empty: {output}")
    block = gemmi.cif.read_file(str(cif)).sole_block()
    sites = block.get_mmcif_category("_atom_site.")
    entities = block.get_mmcif_category("_entity.")
    if not sites or not entities:
        raise ValueError("mmCIF must contain atom_site and entity categories")
    for field in ("id", "label_atom_id", "label_comp_id", "label_asym_id",
                  "label_entity_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z"):
        if field not in sites:
            raise ValueError(f"mmCIF atom_site is missing {field}")
    entity_types = dict(zip(entities["id"], entities["type"]))
    entry_id = block.name.upper()
    grouped = defaultdict(list)
    for i in range(len(sites["id"])):
        entity = _column(sites, "label_entity_id", i)
        if entity not in entity_types:
            raise ValueError(f"Atom references unknown entity {entity!r}")
        if entity_types[entity].lower() != "non-polymer":
            continue
        comp = _column(sites, "label_comp_id", i)
        if components is not None and comp.upper() not in components:
            continue
        key = (entry_id, _column(sites, "pdbx_PDB_model_num", i, "1"),
               _column(sites, "label_asym_id", i), _column(sites, "auth_asym_id", i),
               _column(sites, "auth_seq_id", i), _column(sites, "pdbx_PDB_ins_code", i), comp)
        grouped[key].append({
            "atom_id": _column(sites, "label_atom_id", i), "alt": _column(sites, "label_alt_id", i),
            "element": _column(sites, "type_symbol", i), "x": _column(sites, "Cartn_x", i),
            "y": _column(sites, "Cartn_y", i), "z": _column(sites, "Cartn_z", i),
        })
    cache = {}
    seen = {}
    records = []
    client = client or PDBClient()
    for key, atoms in grouped.items():
        alts = sorted({atom["alt"] for atom in atoms if atom["alt"]}) or [""]
        for alt in alts:
            record = dict(zip(FIELDS[:7], key))
            record.update(alt_id=alt, file="", status="error", detail="")
            try:
                component_id = key[-1].upper()
                if component_id not in cache:
                    ccd_path = (Path(ccd_dir) if ccd_dir else output / "ccd") / f"{component_id}.cif"
                    if ccd_dir is None and not ccd_path.exists():
                        ccd_path = client.download(component_id, "component", "cif", output / "ccd")
                    cache[component_id] = _component(ccd_path)
                selected = [atom for atom in atoms if not atom["alt"] or atom["alt"] == alt]
                molecule, missing, unknown = _molecule(selected, *cache[component_id],
                                                       allow_unknown_charge=allow_unknown_charge)
                notes = []
                if missing:
                    notes.append(f"Missing CCD heavy atoms: {', '.join(missing)}")
                if unknown:
                    notes.append(f"Unknown CCD formal charges assumed zero: {', '.join(unknown)}")
                if pattern is not None and not molecule.HasSubstructMatch(pattern):
                    record.update(status="filtered", detail="; ".join(["SMARTS did not match", *notes]))
                else:
                    signature = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False)
                    if deduplicate and signature in seen:
                        record.update(file=seen[signature], status="duplicate",
                                      detail="; ".join(["Same canonical connectivity as earlier occurrence", *notes]))
                        records.append(record)
                        continue
                    filename = f"ligand_{len(records) + 1:06}.sdf"
                    molecule.SetProp("_Name", f"{entry_id} {component_id}")
                    writer = Chem.SDWriter(str(output / filename))
                    try:
                        writer.write(molecule)
                    finally:
                        writer.close()
                    record.update(file=filename, status="written", detail="; ".join(notes))
                    seen[signature] = filename
            except (OSError, ValueError, RuntimeError, requests.RequestException) as exc:
                record["detail"] = str(exc)
            records.append(record)
    with (output / "manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(records)
    failures = sum(row["status"] == "error" for row in records)
    if failures:
        raise RuntimeError(f"{failures} ligand occurrence(s) failed; see {output / 'manifest.csv'}")
    return records
