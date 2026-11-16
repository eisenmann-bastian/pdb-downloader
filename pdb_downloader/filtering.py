"""Filter SDF molecules with a SMARTS pattern."""

from pathlib import Path

from rdkit import Chem


def _read_matching(source: Path, smarts: str | None) -> list:
    pattern = Chem.MolFromSmarts(smarts) if smarts is not None else None
    if smarts is not None and pattern is None:
        raise ValueError(f"Invalid SMARTS pattern: {smarts}")
    supplier = Chem.SDMolSupplier(str(source), removeHs=False)
    molecules = []
    for i, mol in enumerate(supplier, 1):
        if mol is None:
            raise ValueError(f"Invalid SDF molecule at record {i} in {source}")
        if pattern is None or mol.HasSubstructMatch(pattern):
            molecules.append((i, mol))
    return molecules


def filter_sdf(source: Path, output: Path, smarts: str) -> int:
    source, output = Path(source), Path(output)
    if output.exists():
        raise FileExistsError(f"Output already exists: {output}")
    molecules = _read_matching(source, smarts)
    output.parent.mkdir(parents=True, exist_ok=True)
    writer = Chem.SDWriter(str(output))
    try:
        for _, mol in molecules:
            writer.write(mol)
    finally:
        writer.close()
    return len(molecules)


def split_sdf(source: Path, output: Path, smarts: str | None = None) -> list[Path]:
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(f"Split output directory must be empty: {output}")
    molecules = _read_matching(Path(source), smarts)
    output.mkdir(parents=True, exist_ok=True)
    paths = []
    for index, mol in molecules:
        path = output / f"molecule_{index:06}.sdf"
        writer = Chem.SDWriter(str(path))
        try:
            writer.write(mol)
        finally:
            writer.close()
        paths.append(path)
    return paths
