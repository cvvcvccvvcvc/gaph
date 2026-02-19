"""Main pipeline entry point.

Reads config, sets up logging, loops over genes, calls run_gene.

Usage:
    python pipeline/pipeline.py [config.json]
"""

import json
import sys
from datetime import datetime
from pathlib import Path

from Bio import Entrez
from loguru import logger

import config
from run_gene import run_gene


def init_logging(log_path):
    """Configure loguru handlers for this run."""
    logger.remove()
    logger.add(
        sys.stdout,
        format="<green>{time:HH:mm:ss}</green> | <level>{level: <8}</level> | <level>{message}</level>",
        level="DEBUG",
        colorize=True,
    )
    logger.add(
        str(log_path),
        format="{time:YYYY-MM-DD HH:mm:ss} | {level: <8} | {message}",
        level="DEBUG",
        rotation="10 MB",
        retention="7 days",
    )
    logger.info(f"Logging initialized. Log file: {log_path}")


def main():
    # Load config
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    with open(config_path) as f:
        cfg = json.load(f)
    config.init(cfg)

    # Set Entrez credentials
    Entrez.email = config.ENTREZ_EMAIL or None
    Entrez.api_key = config.ENTREZ_API_KEY or None
    if not Entrez.email:
        logger.warning("ENTREZ_EMAIL is not set; NCBI may warn or throttle anonymous requests")
    if not Entrez.api_key:
        logger.warning("ENTREZ_API_KEY is not set; using unauthenticated NCBI rate limits")

    # Create run directory
    run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_dir = config.RUNS_DIR / f"run_{run_id}"
    run_dir.mkdir(parents=True, exist_ok=True)

    init_logging(run_dir / "pipeline.log")

    # Save run params
    run_params = {
        "run_id": run_id,
        "run_dir": str(run_dir),
        "gene_ids": cfg["gene_ids"],
        "hitlist_size": cfg.get("hitlist_size", 5000),
        "config_path": str(config_path),
        "started_at": datetime.now().isoformat(timespec="seconds"),
    }
    with open(run_dir / "run_params.json", "w") as f:
        json.dump(run_params, f, indent=2)

    logger.info(f"Run directory: {run_dir}")
    logger.info(f"Gene IDs: {cfg['gene_ids']}")

    # Process each gene (ortholog source: NCBI Datasets -> BLAST fallback)
    for gene_id in cfg["gene_ids"]:
        gene_dir = run_dir / f"gene_{gene_id}"
        try:
            run_gene(gene_id, gene_dir, cfg)
        except Exception:
            logger.exception(f"Failed to process gene {gene_id}")
            continue

    logger.success("Pipeline run complete")


if __name__ == "__main__":
    main()
