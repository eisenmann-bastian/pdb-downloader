# Project instructions

## Structure and conventions

- Keep reusable logic in `pdb_downloader/` and command-line handling in `cli.py`.
- Configure `molecular_utilities.logging` once at the entry point; use child
  loggers under the project namespace in other modules. Do not configure
  logging as an import side effect.
- Keep CLI, logging namespace, tests and documentation consistent.

## Changes and validation

- Preserve existing behavior unless a change is requested. Report errors
  clearly; do not silently ignore failures.
- Add or update tests for behavior changes and keep the README and examples
  consistent with the public CLI and API.
- Run `python -m unittest discover -s tests -v` in an installed environment.
  Run `./examples/run-example.sh` when changing the CLI or setup workflow.
- Do not overwrite an existing virtual environment or commit secrets, generated
  output or local environment files.
