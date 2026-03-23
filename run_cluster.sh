#!/usr/bin/env bash
#SBATCH --job-name=gaph
#SBATCH --partition=main
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/mnt/tank/scratch/%u/GAPH/runs/slurm-%j.out
#SBATCH --error=/mnt/tank/scratch/%u/GAPH/runs/slurm-%j.err

set -euo pipefail
umask 077

cd /nfs/home/$USER/repos/gaph
export PATH="$HOME/bin:$PATH"
export MAMBA_ROOT_PREFIX="$HOME/micromamba"

: "${ENTREZ_EMAIL:?ENTREZ_EMAIL is required}"
: "${ENTREZ_API_KEY:?ENTREZ_API_KEY is required}"

SECRETS_ENV="$(mktemp /tmp/gaph_secrets.XXXXXX)"
RUNTIME_CFG="$(mktemp /tmp/gaph_cfg.XXXXXX.json)"
cleanup() { rm -f "$SECRETS_ENV" "$RUNTIME_CFG"; }
trap cleanup EXIT

cat > "$SECRETS_ENV" <<EOF
ENTREZ_EMAIL=$ENTREZ_EMAIL
ENTREZ_API_KEY=$ENTREZ_API_KEY
EOF

CONFIG_SRC="${1:-config.json}"
python3 - "$SECRETS_ENV" "$CONFIG_SRC" "$RUNTIME_CFG" <<'PY'
import json, sys
from pathlib import Path
env_file, src, dst = sys.argv[1], sys.argv[2], sys.argv[3]
with open(src) as f:
    cfg = json.load(f)
cfg["env_file"] = env_file
cfg["source_config_path"] = src
config_stem = Path(src).stem
if not cfg.get("run_name"):
    cfg["run_name"] = config_stem
with open(dst, "w") as f:
    json.dump(cfg, f, indent=2)
    f.write("\n")
PY

bash run.sh "$RUNTIME_CFG"
