# PDB Downloader

Independent CLI commands and a Python API for searching RCSB PDB, downloading
entries and chemical-component definitions, extracting bound non-polymer
occurrences, and filtering or splitting SDF files. Nextflow and PyMOL are not
required.

## Installation

Requires Python 3.10 or later, Git, and network access to install the shared
logging dependency. RDKit and Gemmi are installed as Python dependencies.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -e .
pdb-downloader --help
```

Alternatively use `./examples/setup-venv.sh`. Existing environments are not
overwritten. The commands can also be run with `python -m pdb_downloader`.

## Search

```sh
pdb-downloader search --protein hemoglobin --organism 'Homo sapiens' \
  --max-resolution 2.0 --max-results 20 --output hits.json
pdb-downloader search --smiles 'CCO' --match graph-exact \
  --return-type component --max-results 10
```

Search prints a JSON array of identifiers to stdout, or writes it to a new
file with `--output`. Metadata criteria apply to entries: protein is a
description phrase, organism is an exact scientific name, and resolution is
an upper bound in angstroms. Multiple criteria are combined with AND at the
entry level; they do not necessarily describe the same entity within an
entry. Chemical SMILES search can return entries or chemical-component
definitions. The default is relaxed substructure matching; `--match` exposes
the other supported RCSB graph match modes. SMARTS is used for local
filtering, not as the remote SMILES query.

Search retrieves pages until complete unless `--max-results` is specified.
When more than 10,000 hits are requested it uses RCSB's `return_all_hits`,
which can consume substantial time and memory. Broad queries can time out
at RCSB; errors are reported rather than silently truncating the result set.

## Download

```sh
pdb-downloader download 4HHB 1CRN --kind entry --format cif --output structures
pdb-downloader download ATP --kind component --format sdf --output components
pdb-downloader download ATP --kind component --format cif --output components
pdb-downloader download PRD_000001 --kind bird --format cif --output components
```

Entries support CIF and PDB; chemical components support CIF and **ideal**
SDF, and BIRD `PRD_` definitions support CIF. Chemical searches can return
both component and BIRD IDs; select the appropriate download kind.
An ideal component SDF is not the bound geometry observed in an entry.
Existing output files are not overwritten unless `--overwrite` is given.
For multiple IDs, failed downloads are reported, other IDs are attempted,
and the command exits nonzero if any failed.

For a runnable example, download human deoxyhemoglobin (PDB ID `4HHB`) in
PDB format:

```sh
./examples/run-example-download-protein.sh
```

The script creates or reuses the project `.venv` and prints a unique output
temporary directory containing `4HHB.pdb`. Results are retained under
`TMPDIR` (or `/tmp` if unset) and are not automatically deleted. Set
`EXAMPLE_OUTPUT_DIR` to use another parent directory. Network access to
RCSB PDB is required.

## Extract bound ligands

```sh
pdb-downloader extract structures/4HHB.cif --component HEM \
  --allow-unknown-charge --output extracted
```

Extraction accepts an mmCIF entry and selects **non-polymer entities**, not
water or polymer/branched entities. `--component` may be repeated; without
it all non-polymer components, including ions, are considered. For each
model and residue occurrence, separate alternate conformers include shared
atoms and the chosen alternate atoms. Each result gets its experimental
coordinates and chemical-component dictionary (CCD) bond orders, not
distance-inferred bonds. CCD definitions are downloaded from RCSB into the
output directory, or supplied offline with `--ccd-dir` containing files named
`COMPONENT_ID.cif`.

The output directory must be empty. `manifest.csv` records the entry, model,
label and author chain IDs, author residue number, insertion code, component,
alternate ID, file, status and details for every occurrence. Filtered and
failed occurrences also appear. Extraction continues after per-occurrence
errors and exits nonzero if any fail; inspect the manifest. Existing output
is never overwritten.

```sh
pdb-downloader extract structures/4HHB.cif --component HEM \
  --smarts '[Fe]' --deduplicate --allow-unknown-charge --output unique-hemes
```

By default all occurrences are written. `--smarts` selects molecules locally;
`--deduplicate` writes only the first molecule for each canonical connectivity
and retains a manifest row pointing to that file for every duplicate. It
ignores geometry and stereochemistry, so use it only if those differences do
not matter.

**Chemical limitations:** the CCD occasionally leaves formal charge unknown
(including the Fe atom in HEM). By default those occurrences fail.
`--allow-unknown-charge` explicitly assumes zero and lists affected atoms in
the manifest. Missing observed heavy atoms are listed there as well; the
result may not represent the complete component. Missing hydrogens are not
reconstructed from experimental coordinates, and stereochemistry is not
assigned from CCD or coordinates. Bonds to other residues, metal
coordination, covalent links, branched ligands and biological-assembly
transformations are not reconstructed. SMARTS matching therefore acts on
the reconstructed individual component, without stereochemical matching.
Review the structure and manifest before using outputs in chemical or
docking calculations.

## Filter or split SDF

```sh
pdb-downloader filter input.sdf --smarts '[#6]-[F,Cl,Br,I]' --output matches.sdf
pdb-downloader split input.sdf --smarts '[#6]' --output split-molecules
```

Filtering writes one SDF containing matched molecules; splitting writes one
file per matching input record, with its original record number in the file
name. `--smarts` is optional for `split`. Invalid SDF records or SMARTS are
errors. Filtering does not overwrite an existing file, and the split output
directory must be empty.

## Python API and development

```python
from pathlib import Path
from pdb_downloader import PDBClient, chemical_query, extract_ligands

client = PDBClient()
ids = client.search(chemical_query("CCO", "graph-exact"), max_results=10)
client.download("4HHB", "entry", "cif", Path("structures"))
records = extract_ligands(Path("structures/4HHB.cif"), Path("extracted"),
                          components={"HEM"}, allow_unknown_charge=True)
```

`metadata_query`, `filter_sdf`, and `split_sdf` are also public. Network
errors and invalid inputs raise exceptions for API callers. CLI diagnostic
logging uses `molecular-utilities` on stderr.

Run `python -m unittest discover -s tests -v`. `./examples/run-example.sh`
performs a small live metadata search. The older `pdb-query-ligand/` and
`pdb-query-receptor/` directories are retained temporarily for migration
reference; they are not dependencies of the new commands.

## License

Distributed under the [BSD 3-Clause License](LICENSE). RCSB PDB data and
third-party dependencies are subject to their own terms.
