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

output_root=${EXAMPLE_OUTPUT_DIR:-${TMPDIR:-/tmp}}
mkdir -p "$output_root"
output_root=$(CDPATH= cd -- "$output_root" && pwd)
output=$(mktemp -d "$output_root/hemoglobin.XXXXXXXX")
printf 'Output directory: %s\n' "$output"
"$venv/bin/python" -m pdb_downloader download 4HHB --kind entry --format pdb --output "$output"
