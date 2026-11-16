"""Command-line tools for querying and downloading RCSB PDB data."""

import argparse
import json
from pathlib import Path

import requests
from molecular_utilities.logging import configure_logging, get_logger

from .client import PDBClient, chemical_query, metadata_query
from .extract import extract_ligands
from .filtering import filter_sdf, split_sdf

logger = get_logger("cli", namespace="pdb_downloader")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)

    search = commands.add_parser("search", help="Search chemical structures and/or PDB metadata")
    search.add_argument("--smiles")
    search.add_argument("--match", default="sub-struct-graph-relaxed",
                        choices=["graph-exact", "graph-relaxed",
                                 "sub-struct-graph-exact", "sub-struct-graph-relaxed"])
    search.add_argument("--protein", help="Protein description phrase")
    search.add_argument("--organism", help="Scientific organism name")
    search.add_argument("--max-resolution", type=float, help="Maximum resolution in angstroms")
    search.add_argument("--return-type", choices=["entry", "component"], default="entry")
    search.add_argument("--page-size", type=int, default=1000)
    search.add_argument("--max-results", type=int)
    search.add_argument("--output", type=Path, help="JSON file; stdout when omitted")

    download = commands.add_parser("download", help="Download entries, components or BIRD definitions by ID")
    download.add_argument("ids", nargs="+")
    download.add_argument("--kind", choices=["entry", "component", "bird"], required=True)
    download.add_argument("--format", choices=["cif", "pdb", "sdf"], required=True)
    download.add_argument("--output", type=Path, required=True)
    download.add_argument("--overwrite", action="store_true")

    extract = commands.add_parser("extract", help="Extract bound non-polymer ligands from mmCIF")
    extract.add_argument("cif", type=Path)
    extract.add_argument("--output", type=Path, required=True)
    extract.add_argument("--component", action="append", help="Include component ID; may be repeated")
    extract.add_argument("--smarts", help="Keep occurrences matching this SMARTS pattern")
    extract.add_argument("--ccd-dir", type=Path, help="Directory containing component CIF files")
    extract.add_argument("--deduplicate", action="store_true",
                         help="Write one SDF per canonical connectivity, preserving all manifest rows")
    extract.add_argument("--allow-unknown-charge", action="store_true",
                         help="Assume zero for unspecified CCD charges and record this in the manifest")

    filt = commands.add_parser("filter", help="Filter SDF records by SMARTS")
    filt.add_argument("sdf", type=Path)
    filt.add_argument("--smarts", required=True)
    filt.add_argument("--output", type=Path, required=True)

    split = commands.add_parser("split", help="Split SDF records into separate files")
    split.add_argument("sdf", type=Path)
    split.add_argument("--smarts", help="Only write records matching this pattern")
    split.add_argument("--output", type=Path, required=True)

    args = parser.parse_args(argv)
    configure_logging(namespace="pdb_downloader")
    try:
        if args.command == "search":
            nodes = []
            if args.smiles is not None:
                nodes.append(chemical_query(args.smiles, args.match))
            if any(value is not None for value in (args.protein, args.organism, args.max_resolution)):
                if args.return_type == "component":
                    raise ValueError("Metadata criteria require --return-type entry")
                nodes.append(metadata_query(args.protein, args.organism, args.max_resolution))
            if not nodes:
                raise ValueError("Provide at least one search criterion")
            query = nodes[0] if len(nodes) == 1 else {
                "type": "group", "logical_operator": "and", "nodes": nodes,
            }
            ids = PDBClient().search(query, "entry" if args.return_type == "entry" else "mol_definition",
                                     page_size=args.page_size, max_results=args.max_results)
            text = json.dumps(ids, indent=2) + "\n"
            if args.output:
                with args.output.open("x", encoding="utf-8") as handle:
                    handle.write(text)
                logger.info("Wrote %d identifiers to %s", len(ids), args.output)
            else:
                print(text, end="")
        elif args.command == "download":
            client = PDBClient()
            failures = []
            for identifier in args.ids:
                try:
                    print(client.download(identifier, args.kind, args.format, args.output, overwrite=args.overwrite))
                except (OSError, ValueError, requests.RequestException) as exc:
                    failures.append(identifier)
                    logger.error("%s: %s", identifier, exc)
            if failures:
                return 1
        elif args.command == "extract":
            records = extract_ligands(args.cif, args.output,
                                      components={x.upper() for x in args.component} if args.component else None,
                                      smarts=args.smarts, ccd_dir=args.ccd_dir,
                                      deduplicate=args.deduplicate,
                                      allow_unknown_charge=args.allow_unknown_charge)
            logger.info("Wrote %d ligand occurrences to %s", sum(r["status"] == "written" for r in records), args.output)
            print(args.output / "manifest.csv")
        elif args.command == "filter":
            count = filter_sdf(args.sdf, args.output, args.smarts)
            logger.info("Wrote %d matching molecules to %s", count, args.output)
            print(args.output)
        else:
            paths = split_sdf(args.sdf, args.output, args.smarts)
            logger.info("Wrote %d molecules to %s", len(paths), args.output)
            for path in paths:
                print(path)
        return 0
    except (OSError, ValueError, RuntimeError, requests.RequestException, KeyError) as exc:
        logger.error("%s", exc)
        return 1
