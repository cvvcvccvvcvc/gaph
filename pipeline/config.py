"""Central configuration for the bioinformatics pipeline.

All settings come from config.json passed at startup.
Secrets (ENTREZ_EMAIL, ENTREZ_API_KEY) live in a separate .env file
whose path is specified by the "env_file" key in config.json.
"""

import shutil
import subprocess
from pathlib import Path

from loguru import logger

REPO_ROOT = Path(__file__).parent.parent.resolve()

# Module globals — set by init()
DATA_DIR = None
RUNS_DIR = None
GNOMAD_CACHE_DIR = None
CLINVAR_VCF = None
CONDA_ENV_NAME = None
CONDA_EXE = None
ENTREZ_EMAIL = ""
ENTREZ_API_KEY = ""


def _resolve_path(p):
    """Resolve path: expand ~, make relative paths relative to REPO_ROOT."""
    p = Path(p).expanduser()
    return p if p.is_absolute() else REPO_ROOT / p


def _load_env_file(path):
    """Parse a .env file and return dict of key=value pairs.

    Supports lines like:
        KEY=VALUE
        export KEY=VALUE
        KEY="VALUE"
    Ignores comments (#) and blank lines.
    """
    env = {}
    path = Path(path)
    if not path.exists():
        logger.warning(f"env_file not found: {path}")
        return env
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("export "):
            line = line[len("export "):]
        if "=" not in line:
            continue
        key, _, value = line.partition("=")
        key = key.strip()
        value = value.strip().strip("\"'")
        env[key] = value
    return env


def _find_conda():
    for exe in ["micromamba", "mamba", "conda"]:
        if shutil.which(exe):
            return exe
    # Check common locations not in PATH
    home = Path.home()
    for exe in ["micromamba", "mamba", "conda"]:
        for loc in [home / "bin" / exe, home / ".local" / "bin" / exe]:
            if loc.exists():
                return str(loc)
    raise RuntimeError("No conda/micromamba/mamba found in PATH or ~/bin")


def init(cfg: dict):
    """Initialize config from parsed JSON dict."""
    global DATA_DIR, RUNS_DIR, GNOMAD_CACHE_DIR, CLINVAR_VCF
    global CONDA_ENV_NAME, CONDA_EXE, ENTREZ_EMAIL, ENTREZ_API_KEY

    # Reset credentials on every init() to avoid stale values across re-initialization.
    ENTREZ_EMAIL = ""
    ENTREZ_API_KEY = ""

    # Load .env file if specified
    env_file = cfg.get("env_file")
    if env_file:
        env_vars = _load_env_file(_resolve_path(env_file))
        ENTREZ_EMAIL = env_vars.get("ENTREZ_EMAIL", "")
        ENTREZ_API_KEY = env_vars.get("ENTREZ_API_KEY", "")

    # Paths
    DATA_DIR = _resolve_path(cfg.get("data_dir", "data"))
    RUNS_DIR = _resolve_path(cfg.get("runs_dir", "runs"))
    GNOMAD_CACHE_DIR = _resolve_path(cfg.get("gnomad_dir", "gnomad_variants"))
    CLINVAR_VCF = DATA_DIR / "clinvar.vcf.gz"

    # Conda
    CONDA_ENV_NAME = cfg.get("conda_env", "bio")
    CONDA_EXE = _find_conda()


def run_bio(cmd):
    """Run a command inside the conda bio environment."""
    logger.debug(f"Running: {cmd}")
    full = f"{CONDA_EXE} run -n {CONDA_ENV_NAME} {cmd}"
    try:
        result = subprocess.run(full, shell=True, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        if e.stdout:
            logger.error(e.stdout.strip())
        if e.stderr:
            logger.error(e.stderr.strip())
        logger.error(f"Command failed (exit={e.returncode}): {cmd}")
        raise
    if result.stdout:
        logger.debug(result.stdout.strip())
    if result.stderr:
        logger.debug(result.stderr.strip())
    return result
