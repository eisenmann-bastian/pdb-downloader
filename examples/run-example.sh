#!/bin/sh
set -eu

root=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
if [ "$#" -ne 0 ]; then
    echo "Usage: $0" >&2
    exit 2
fi
venv="$root/.venv"
if [ ! -e "$venv" ] && [ ! -L "$venv" ]; then
    "$root/examples/setup-venv.sh" "$venv"
fi
if [ ! -f "$venv/bin/activate" ] || [ ! -x "$venv/bin/python" ]; then
    echo "Invalid virtual environment at $venv; inspect it before retrying." >&2
    exit 1
fi
. "$venv/bin/activate"
"$venv/bin/python" -m pdb_downloader search --protein hemoglobin --max-results 5
