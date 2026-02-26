#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${1:-config.json}"

# Ensure ~/bin is in PATH (micromamba install location)
export PATH="$HOME/bin:$PATH"

if command -v python3 &>/dev/null; then
    JSON_PYTHON="python3"
elif command -v python &>/dev/null; then
    JSON_PYTHON="python"
else
    echo "ERROR: python3/python not found to parse config file: $CONFIG" >&2
    exit 1
fi

CONDA_ENV="$("$JSON_PYTHON" - "$CONFIG" <<'PY'
import json
import sys
from pathlib import Path

cfg_path = Path(sys.argv[1]).expanduser()
with cfg_path.open() as f:
    cfg = json.load(f)
print(cfg.get("conda_env", "bio"))
PY
)"

# Resolve conda-compatible runner.
# Using `run -n` avoids accidental use of an active `.venv` python from parent shell.
if command -v micromamba &>/dev/null; then
    CONDA_RUNNER="micromamba"
elif command -v mamba &>/dev/null; then
    CONDA_RUNNER="mamba"
elif command -v conda &>/dev/null; then
    CONDA_RUNNER="conda"
else
    echo "ERROR: No conda/micromamba/mamba found" >&2
    exit 1
fi

export PYTHONPATH="$SCRIPT_DIR/pipeline${PYTHONPATH:+:$PYTHONPATH}"

if "$CONDA_RUNNER" run -n "$CONDA_ENV" python -V &>/dev/null; then
    PIPELINE_PYTHON="python"
elif "$CONDA_RUNNER" run -n "$CONDA_ENV" python3 -V &>/dev/null; then
    PIPELINE_PYTHON="python3"
else
    echo "ERROR: Neither python nor python3 found in env: $CONDA_ENV" >&2
    exit 1
fi

"$CONDA_RUNNER" run -n "$CONDA_ENV" "$PIPELINE_PYTHON" "$SCRIPT_DIR/pipeline/pipeline.py" "$CONFIG"
