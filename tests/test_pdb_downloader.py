import csv
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import requests
from rdkit import Chem

from pdb_downloader.client import PDBClient, chemical_query, metadata_query
from pdb_downloader.cli import main
from pdb_downloader.extract import extract_ligands
from pdb_downloader.filtering import filter_sdf, split_sdf


STRUCTURE = """data_TEST
loop_
_entity.id
_entity.type
1 non-polymer
2 water
loop_
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.auth_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.pdbx_PDB_model_num
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
1 C C1 . LIG A 1 X 5 ? 1 0 0 0
2 O O1 A LIG A 1 X 5 ? 1 1 0 0
3 O O1 B LIG A 1 X 5 ? 1 2 0 0
4 C C1 . LIG B 1 Y 6 ? 1 3 0 0
5 O O1 . LIG B 1 Y 6 ? 1 4 0 0
6 C C1 . LIG A 1 X 5 ? 2 5 0 0
7 O O1 . LIG A 1 X 5 ? 2 6 0 0
8 O O . HOH C 2 Z 7 ? 1 0 0 0
"""
CCD = """data_LIG
loop_
_chem_comp_atom.atom_id
_chem_comp_atom.type_symbol
_chem_comp_atom.charge
C1 C 0
O1 O 0
loop_
_chem_comp_bond.atom_id_1
_chem_comp_bond.atom_id_2
_chem_comp_bond.value_order
C1 O1 sing
"""


def response(data=None, status=200, content=b"data"):
    result = Mock()
    result.status_code = status
    result.json.return_value = data
    result.content = content
    return result


class ClientTests(unittest.TestCase):
    def test_pagination_and_result_limit(self):
        session = Mock()
        session.post.side_effect = [
            response({"total_count": 3, "result_set": [{"identifier": "A"}, {"identifier": "B"}]}),
            response({"total_count": 3, "result_set": [{"identifier": "C"}]}),
        ]
        result = PDBClient(session).search(chemical_query("CCO"), "mol_definition", page_size=2)
        self.assertEqual(result, ["A", "B", "C"])
        self.assertEqual(session.post.call_args_list[1].kwargs["json"]["request_options"]["paginate"],
                         {"start": 2, "rows": 2})
        session.post.side_effect = [response({"total_count": 3, "result_set": [{"identifier": "A"}]})]
        self.assertEqual(PDBClient(session).search(chemical_query("CCO"), max_results=1), ["A"])

    def test_empty_and_incomplete_search(self):
        session = Mock()
        session.post.return_value = response(status=204)
        self.assertEqual(PDBClient(session).search(chemical_query("C")), [])
        session.post.return_value = response({"total_count": 1, "result_set": []})
        with self.assertRaisesRegex(ValueError, "empty page"):
            PDBClient(session).search(chemical_query("C"))
        with self.assertRaises(ValueError):
            chemical_query("C", "exact")
        with self.assertRaises(ValueError):
            metadata_query(resolution=float("nan"))

    def test_large_result_set_uses_return_all_hits(self):
        session = Mock()
        all_results = [{"identifier": str(i)} for i in range(10001)]
        session.post.side_effect = [
            response({"total_count": 10001, "result_set": all_results[:2]}),
            response({"total_count": 10001, "result_set": all_results}),
        ]
        result = PDBClient(session).search(chemical_query("C"), page_size=2)
        self.assertEqual(len(result), 10001)
        self.assertEqual(session.post.call_args_list[1].kwargs["json"]["request_options"],
                         {"return_all_hits": True})

    def test_metadata_combination(self):
        query = metadata_query("hemoglobin", "Homo sapiens", 2.0)
        self.assertEqual(query["logical_operator"], "and")
        self.assertEqual(len(query["nodes"]), 3)
        self.assertEqual(query["nodes"][2]["parameters"]["value"], 2.0)

    def test_download_url_safety_and_overwrite(self):
        session = Mock()
        session.get.return_value = response(content=b"molecule")
        with tempfile.TemporaryDirectory() as tmp:
            path = PDBClient(session).download("atp", "component", "sdf", Path(tmp))
            self.assertEqual(path.read_bytes(), b"molecule")
            self.assertIn("ATP_ideal.sdf", session.get.call_args.args[0])
            with self.assertRaises(FileExistsError):
                PDBClient(session).download("atp", "component", "sdf", Path(tmp))
            with self.assertRaises(ValueError):
                PDBClient(session).download("../x", "entry", "cif", Path(tmp))
            with self.assertRaises(ValueError):
                PDBClient(session).download("atp", "entry", "sdf", Path(tmp))
            bird = PDBClient(session).download("PRD_000001", "bird", "cif", Path(tmp))
            self.assertEqual(bird.name, "PRD_000001.cif")
            self.assertIn("/birds/download/PRD_000001.cif", session.get.call_args.args[0])


class ExtractTests(unittest.TestCase):
    def setUp(self):
        temp = tempfile.TemporaryDirectory()
        self.addCleanup(temp.cleanup)
        self.root = Path(temp.name)
        self.structure = self.root / "structure.cif"
        self.structure.write_text(STRUCTURE)
        self.ccd = self.root / "ccd"
        self.ccd.mkdir()
        (self.ccd / "LIG.cif").write_text(CCD)

    def test_occurrences_alternates_and_coordinates(self):
        out = self.root / "out"
        records = extract_ligands(self.structure, out, ccd_dir=self.ccd)
        self.assertEqual(len(records), 4)
        self.assertEqual([r["alt_id"] for r in records], ["A", "B", "", ""])
        self.assertEqual([r["model"] for r in records], ["1", "1", "1", "2"])
        self.assertEqual([r["auth_asym_id"] for r in records], ["X", "X", "Y", "X"])
        with (out / "manifest.csv").open() as handle:
            self.assertEqual(len(list(csv.DictReader(handle))), 4)
        for i, record in enumerate(records):
            mol = Chem.SDMolSupplier(str(out / record["file"]))[0]
            self.assertEqual(mol.GetNumAtoms(), 2)
            self.assertEqual(mol.GetNumBonds(), 1)
            self.assertAlmostEqual(mol.GetConformer().GetAtomPosition(1).x, [1, 2, 4, 6][i])

    def test_deduplicate_keeps_provenance(self):
        out = self.root / "out"
        records = extract_ligands(self.structure, out, ccd_dir=self.ccd, deduplicate=True)
        self.assertEqual([r["status"] for r in records], ["written", "duplicate", "duplicate", "duplicate"])
        self.assertEqual(len({r["file"] for r in records}), 1)
        self.assertEqual(len(list(out.glob("*.sdf"))), 1)

    def test_filter_and_failed_component_are_reported(self):
        out = self.root / "filtered"
        records = extract_ligands(self.structure, out, ccd_dir=self.ccd, smarts="[N]")
        self.assertTrue(all(r["status"] == "filtered" for r in records))
        out = self.root / "failed"
        with self.assertRaisesRegex(RuntimeError, "4 ligand"):
            extract_ligands(self.structure, out, ccd_dir=self.root / "missing")
        with (out / "manifest.csv").open() as handle:
            self.assertTrue(all(r["status"] == "error" for r in csv.DictReader(handle)))
        with self.assertRaises(ValueError):
            extract_ligands(self.structure, self.root / "bad", smarts="[")
        with self.assertRaises(FileExistsError):
            extract_ligands(self.structure, out, ccd_dir=self.ccd)

    def test_remote_failure_still_writes_manifest(self):
        client = Mock()
        client.download.side_effect = requests.ConnectionError("offline")
        with self.assertRaisesRegex(RuntimeError, "failed"):
            extract_ligands(self.structure, self.root / "offline", client=client)
        self.assertIn("offline", (self.root / "offline/manifest.csv").read_text())

    def test_unknown_charge_requires_explicit_opt_in(self):
        (self.ccd / "LIG.cif").write_text(CCD.replace("O1 O 0", "O1 O ?"))
        with self.assertRaises(RuntimeError):
            extract_ligands(self.structure, self.root / "strict", ccd_dir=self.ccd)
        records = extract_ligands(self.structure, self.root / "allowed", ccd_dir=self.ccd,
                                  allow_unknown_charge=True)
        self.assertTrue(all("assumed zero: O1" in r["detail"] for r in records))

    def test_missing_bonds_are_reported(self):
        (self.ccd / "LIG.cif").write_text(CCD.replace("C1 O1 sing\n", ""))
        with self.assertRaises(RuntimeError):
            extract_ligands(self.structure, self.root / "no-bonds", ccd_dir=self.ccd)
        self.assertIn("No CCD bonds connect", (self.root / "no-bonds/manifest.csv").read_text())


class FilterAndCliTests(unittest.TestCase):
    def test_filter_preserves_matched_sdf_and_invalid_records_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            src = root / "input.sdf"
            writer = Chem.SDWriter(str(src))
            writer.write(Chem.MolFromSmiles("CCO"))
            writer.write(Chem.MolFromSmiles("CC"))
            writer.close()
            out = root / "output.sdf"
            self.assertEqual(filter_sdf(src, out, "[O]"), 1)
            self.assertEqual(len(Chem.SDMolSupplier(str(out))), 1)
            split = split_sdf(src, root / "split", "[C]")
            self.assertEqual([x.name for x in split], ["molecule_000001.sdf", "molecule_000002.sdf"])
            with self.assertRaises(FileExistsError):
                split_sdf(src, root / "split")
            with self.assertRaises(FileExistsError):
                filter_sdf(src, out, "[O]")
            with self.assertRaises(ValueError):
                filter_sdf(src, root / "bad.sdf", "[")
            src.write_text("invalid\n$$$$\n")
            with self.assertRaises(ValueError):
                filter_sdf(src, root / "invalid.sdf", "[O]")

    def test_cli_search_and_download_failures(self):
        with patch("pdb_downloader.cli.PDBClient") as client:
            client.return_value.search.return_value = ["4HHB"]
            with tempfile.TemporaryDirectory() as tmp:
                output = Path(tmp) / "hits.json"
                self.assertEqual(main(["search", "--protein", "hemoglobin", "--output", str(output)]), 0)
                self.assertEqual(json.loads(output.read_text()), ["4HHB"])
                self.assertEqual(main(["search", "--protein", "x", "--return-type", "component"]), 1)
                client.return_value.download.side_effect = requests.ConnectionError("offline")
                self.assertEqual(main(["download", "4HHB", "--kind", "entry", "--format", "cif",
                                       "--output", tmp]), 1)


if __name__ == "__main__":
    unittest.main()
