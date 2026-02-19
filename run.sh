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

# Activate conda environment
if command -v micromamba &>/dev/null; then
    eval "$(micromamba shell hook --shell bash)"
    micromamba activate "$CONDA_ENV"
elif command -v mamba &>/dev/null; then
    eval "$(mamba shell hook --shell bash)"
    mamba activate "$CONDA_ENV"
elif command -v conda &>/dev/null; then
    eval "$(conda shell.bash hook)"
    conda activate "$CONDA_ENV"
else
    echo "ERROR: No conda/micromamba/mamba found" >&2
    exit 1
fi

export PYTHONPATH="$SCRIPT_DIR/pipeline${PYTHONPATH:+:$PYTHONPATH}"

if command -v python &>/dev/null; then
    PIPELINE_PYTHON="python"
elif command -v python3 &>/dev/null; then
    PIPELINE_PYTHON="python3"
else
    echo "ERROR: No python/python3 found after activating conda env: $CONDA_ENV" >&2
    exit 1
fi

"$PIPELINE_PYTHON" "$SCRIPT_DIR/pipeline/pipeline.py" "$CONFIG"
