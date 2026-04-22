"""Main pipeline entry point.

Reads config, sets up logging, loops over genes, calls run_gene.

Usage:
    python pipeline/pipeline.py [config.json]
"""

import json
import re
import shutil
import sys
import traceback
from datetime import datetime
from pathlib import Path

from Bio import Entrez
from loguru import logger

import config
from orthologs import get_source
from run_compaction import (
    compact_run_in_place,
    load_compacted_gene_status,
    run_has_compacted_snapshot,
)
from run_gene import run_gene

GENE_SUCCESS_ARTIFACT = "gene_snps_annotated.vcf"
GENE_FAILURE_ARTIFACT = "failure.json"


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


def _validate_percent(name: str, value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"bam_filtering.{name} must be a number in [0, 100], got: {value!r}")
    value = float(value)
    if value < 0 or value > 100:
        raise ValueError(f"bam_filtering.{name} must be in [0, 100], got: {value}")
    return value


def _validate_positive_int_field(name: str, value) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise TypeError(f"{name} must be a positive integer, got: {value!r}")
    return value


def _validate_positive_number_field(name: str, value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or float(value) <= 0:
        raise TypeError(f"{name} must be a positive number, got: {value!r}")
    return float(value)


def _validate_fraction_field(name: str, value) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError(f"{name} must be a number in [0, 1], got: {value!r}")
    value = float(value)
    if value < 0 or value > 1:
        raise ValueError(f"{name} must be in [0, 1], got: {value}")
    return value


def _validate_phred_field(name: str, value) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{name} must be an integer in [0, 93], got: {value!r}")
    if value < 0 or value > 93:
        raise ValueError(f"{name} must be in [0, 93], got: {value}")
    return value


def validate_bam_filtering_cfg(cfg: dict) -> dict:
    """Validate and normalize the bam_filtering config section."""
    if "bam_filtering" not in cfg:
        raise ValueError(
            "Missing required 'bam_filtering' block in config.json "
            "(filtering is ON by default and must be configured explicitly)"
        )

    raw = cfg["bam_filtering"]
    if not isinstance(raw, dict):
        raise TypeError(f"'bam_filtering' must be an object, got: {type(raw).__name__}")

    enabled = raw.get("enabled", True)
    if not isinstance(enabled, bool):
        raise TypeError(f"bam_filtering.enabled must be boolean, got: {enabled!r}")

    def _read_bool(name: str, default: bool = True) -> bool:
        value = raw.get(name, default)
        if not isinstance(value, bool):
            raise TypeError(f"bam_filtering.{name} must be boolean, got: {value!r}")
        return value

    validated = {
        "enabled": enabled,
        "wrong_strand": _read_bool("wrong_strand", True),
        "lis": _read_bool("lis", True),
        "overlap": _read_bool("overlap", True),
    }

    required_thresholds = [
        "min_mapped_pct_of_generated",
        "max_pct_filtered",
        "min_kept_pct_of_reference",
    ]
    for key in required_thresholds:
        if enabled and key not in raw:
            raise ValueError(
                f"bam_filtering.{key} is required when bam_filtering.enabled=true"
            )

        if key in raw:
            value = raw[key]
            if value is None:
                if enabled:
                    raise ValueError(
                        f"bam_filtering.{key} cannot be null when bam_filtering.enabled=true"
                    )
                validated[key] = None
            else:
                validated[key] = _validate_percent(key, value)
        else:
            validated[key] = None

    return validated


def validate_read_generation_cfg(cfg: dict) -> dict:
    """Validate and normalize read generation parameters."""
    raw = cfg.get("read_generation", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"'read_generation' must be an object, got: {type(raw).__name__}")

    return {
        "pseudo_read_phred": _validate_phred_field(
            "read_generation.pseudo_read_phred",
            raw.get("pseudo_read_phred", 30),
        ),
        "read_len": _validate_positive_int_field(
            "read_generation.read_len",
            raw.get("read_len", 75),
        ),
        "step": _validate_positive_int_field(
            "read_generation.step",
            raw.get("step", 35),
        ),
    }


def validate_variant_calling_cfg(cfg: dict) -> dict:
    """Validate and normalize variant calling parameters."""
    raw = cfg.get("variant_calling", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"'variant_calling' must be an object, got: {type(raw).__name__}")

    return {
        "min_var_freq": _validate_fraction_field(
            "variant_calling.min_var_freq",
            raw.get("min_var_freq", 0.2),
        ),
        "min_coverage": _validate_positive_int_field(
            "variant_calling.min_coverage",
            raw.get("min_coverage", 8),
        ),
        "min_reads2": _validate_positive_int_field(
            "variant_calling.min_reads2",
            raw.get("min_reads2", 2),
        ),
    }


def validate_ortholog_selection_cfg(cfg: dict) -> dict:
    """Validate and normalize ortholog selection parameters."""
    raw = cfg.get("ortholog_selection", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"'ortholog_selection' must be an object, got: {type(raw).__name__}")

    scope = raw.get("scope", "all")
    if not isinstance(scope, str):
        raise TypeError(f"ortholog_selection.scope must be a string, got: {scope!r}")
    scope = scope.strip().lower()
    if not scope:
        scope = "all"

    return {"scope": scope}


def validate_cache_cfg(cfg: dict) -> dict:
    """Validate and normalize optional cache configuration."""
    raw = cfg.get("cache", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"'cache' must be an object, got: {type(raw).__name__}")

    enabled = raw.get("enabled", True)
    if not isinstance(enabled, bool):
        raise TypeError(f"cache.enabled must be boolean, got: {enabled!r}")

    def _read_bool(name: str, default: bool) -> bool:
        value = raw.get(name, default)
        if not isinstance(value, bool):
            raise TypeError(f"cache.{name} must be boolean, got: {value!r}")
        return value

    batch_size = raw.get("ortholog_batch_size", 10)
    if isinstance(batch_size, bool) or not isinstance(batch_size, int) or batch_size <= 0:
        raise TypeError(
            f"cache.ortholog_batch_size must be a positive integer, got: {batch_size!r}"
        )

    validated = {
        "enabled": enabled,
        "orthologs_ncbi": _read_bool("orthologs_ncbi", True),
        "orthologs_blast": _read_bool("orthologs_blast", True),
        "gene_seq": _read_bool("gene_seq", True),
        "gnomad": _read_bool("gnomad", True),
        "ncbi_to_ensembl": _read_bool("ncbi_to_ensembl", True),
        "ortholog_batch_size": batch_size,
    }

    for path_key in [
        "ortholog_cache_dir",
        "gene_seq_cache_dir",
        "ncbi_to_ensembl_cache_file",
    ]:
        if path_key in raw:
            value = raw[path_key]
            if not isinstance(value, str) or not value.strip():
                raise TypeError(f"cache.{path_key} must be a non-empty string path, got: {value!r}")
            validated[path_key] = value

    return validated


def validate_output_compaction_cfg(cfg: dict) -> dict:
    """Validate and normalize optional post-run output compaction."""
    raw = cfg.get("output_compaction", {})
    if raw is None:
        raw = {}
    if not isinstance(raw, dict):
        raise TypeError(f"'output_compaction' must be an object, got: {type(raw).__name__}")

    enabled = raw.get("enabled", False)
    if not isinstance(enabled, bool):
        raise TypeError(f"output_compaction.enabled must be boolean, got: {enabled!r}")

    return {"enabled": enabled}


def _write_json(path: Path, payload) -> None:
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)


def _append_jsonl(path: Path, payload) -> None:
    with open(path, "a") as f:
        f.write(json.dumps(payload, ensure_ascii=False) + "\n")


def _validate_resume_run_dir_cfg(cfg: dict) -> str | None:
    value = cfg.get("resume_run_dir")
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise TypeError(f"resume_run_dir must be a non-empty string path, got: {value!r}")
    return value.strip()


def _validate_run_name_cfg(cfg: dict) -> str | None:
    value = cfg.get("run_name")
    if value is None:
        return None
    if not isinstance(value, str):
        raise TypeError(f"run_name must be a string, got: {value!r}")
    value = value.strip()
    if not value:
        raise ValueError("run_name must not be empty when provided")
    return value


def _resolve_repo_path(path_str: str) -> Path:
    candidate = Path(path_str).expanduser()
    if candidate.is_absolute():
        return candidate.resolve()
    return (config.REPO_ROOT / candidate).resolve()


def _read_run_params_gene_ids(run_dir: Path) -> list[int] | None:
    run_params_path = run_dir / "run_params.json"
    if not run_params_path.exists():
        return None
    with open(run_params_path) as f:
        payload = json.load(f)
    raw_gene_ids = payload.get("gene_ids")
    if raw_gene_ids is None:
        return None
    if not isinstance(raw_gene_ids, list):
        raise TypeError(
            f"{run_params_path} has invalid gene_ids type: expected list, got {type(raw_gene_ids).__name__}"
        )
    return [int(gene_id) for gene_id in raw_gene_ids]


def _gene_status(run_dir: Path, gene_id: int) -> str:
    gene_dir = run_dir / f"gene_{gene_id}"
    success_path = gene_dir / GENE_SUCCESS_ARTIFACT
    failure_path = gene_dir / GENE_FAILURE_ARTIFACT

    if success_path.exists() and success_path.is_file() and success_path.stat().st_size > 0:
        return "success"
    if failure_path.exists() and failure_path.is_file() and failure_path.stat().st_size > 0:
        return "failed"
    if gene_dir.exists():
        return "incomplete"

    compacted_status = load_compacted_gene_status(run_dir, gene_id)
    if compacted_status in {"success", "failed", "incomplete", "missing"}:
        return compacted_status

    return "missing"


def _build_resume_plan(run_dir: Path, gene_ids: list[int]) -> dict:
    counts = {"success": 0, "failed": 0, "incomplete": 0, "missing": 0}
    start_index = None
    start_status = None

    for idx, gene_id in enumerate(gene_ids):
        status = _gene_status(run_dir, gene_id)
        counts[status] += 1
        if start_index is None and status in {"incomplete", "missing"}:
            start_index = idx
            start_status = status

    return {
        "start_index": start_index,
        "start_status": start_status,
        "counts": counts,
    }


def _slugify_run_name(value: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9_-]+", "_", value.strip())
    slug = re.sub(r"_+", "_", slug).strip("_-")
    if not slug:
        raise ValueError(
            f"run_name must contain at least one ASCII letter or digit after sanitization, got: {value!r}"
        )
    return slug


def _derive_run_name(config_path: str, cfg: dict) -> str | None:
    explicit = cfg.get("run_name")
    if explicit:
        return _slugify_run_name(explicit)

    source_config_path = cfg.get("source_config_path") or config_path
    config_stem = Path(source_config_path).stem
    return _slugify_run_name(config_stem)


def _create_new_run_dir(run_name: str | None) -> tuple[str, Path]:
    base_run_id = run_name or "config"

    for suffix in range(1, 10_000):
        run_id = base_run_id if suffix == 1 else f"{base_run_id}__{suffix}"
        run_dir = config.RUNS_DIR / f"run_{run_id}"
        if run_dir.exists():
            continue
        run_dir.mkdir(parents=True, exist_ok=False)
        return run_id, run_dir

    raise RuntimeError(f"Could not allocate unique run directory for base name: {base_run_id}")


def main():
    # Load config
    config_path = sys.argv[1] if len(sys.argv) > 1 else "config.json"
    with open(config_path) as f:
        cfg = json.load(f)
    cfg["hitlist_size"] = _validate_positive_int_field(
        "hitlist_size", cfg.get("hitlist_size", 5000)
    )
    cfg["blast_expect"] = _validate_positive_number_field(
        "blast_expect", cfg.get("blast_expect", 10.0)
    )
    cfg["read_generation"] = validate_read_generation_cfg(cfg)
    cfg["variant_calling"] = validate_variant_calling_cfg(cfg)
    cfg["ortholog_selection"] = validate_ortholog_selection_cfg(cfg)
    cfg["resume_run_dir"] = _validate_resume_run_dir_cfg(cfg)
    cfg["run_name"] = _validate_run_name_cfg(cfg)
    cfg["bam_filtering"] = validate_bam_filtering_cfg(cfg)
    cfg["cache"] = validate_cache_cfg(cfg)
    cfg["output_compaction"] = validate_output_compaction_cfg(cfg)
    config.init(cfg)
    derived_run_name = _derive_run_name(config_path, cfg)

    # Set Entrez credentials
    Entrez.email = config.ENTREZ_EMAIL or None
    Entrez.api_key = config.ENTREZ_API_KEY or None
    if not Entrez.email:
        logger.warning("ENTREZ_EMAIL is not set; NCBI may warn or throttle anonymous requests")
    if not Entrez.api_key:
        logger.warning("ENTREZ_API_KEY is not set; using unauthenticated NCBI rate limits")

    resume_plan = None
    resume_run_dir = None
    resume_mode = False
    start_idx = 0
    if cfg["resume_run_dir"]:
        resume_run_dir = _resolve_repo_path(cfg["resume_run_dir"])
        if not resume_run_dir.exists() or not resume_run_dir.is_dir():
            raise FileNotFoundError(f"resume_run_dir does not exist or is not a directory: {resume_run_dir}")

        existing_gene_ids = _read_run_params_gene_ids(resume_run_dir)
        current_gene_ids = [int(gene_id) for gene_id in cfg["gene_ids"]]
        if existing_gene_ids is not None and existing_gene_ids != current_gene_ids:
            raise ValueError(
                "resume_run_dir gene_ids do not match current config.\n"
                f"run_params gene_ids: {existing_gene_ids}\n"
                f"config gene_ids: {current_gene_ids}"
            )

        resume_plan = _build_resume_plan(resume_run_dir, cfg["gene_ids"])
        if run_has_compacted_snapshot(resume_run_dir) and resume_plan["start_index"] is not None:
            raise ValueError(
                "Cannot resume a compacted run with non-terminal genes. "
                f"Run: {resume_run_dir}; progress snapshot: {resume_plan['counts']}"
            )
        if resume_plan["start_index"] is not None:
            run_dir = resume_run_dir
            run_id = run_dir.name[4:] if run_dir.name.startswith("run_") else run_dir.name
            start_idx = int(resume_plan["start_index"])
            resume_mode = True
        else:
            run_id, run_dir = _create_new_run_dir(derived_run_name)
    else:
        run_id, run_dir = _create_new_run_dir(derived_run_name)

    init_logging(run_dir / "pipeline.log")

    if not resume_mode:
        # Save run params for a newly created run directory.
        run_params = {
            "run_id": run_id,
            "run_dir": str(run_dir),
            "gene_ids": cfg["gene_ids"],
            "hitlist_size": cfg["hitlist_size"],
            "blast_expect": cfg["blast_expect"],
            "read_generation": cfg["read_generation"],
            "variant_calling": cfg["variant_calling"],
            "ortholog_selection": cfg["ortholog_selection"],
            "run_name": derived_run_name,
            "resume_run_dir": str(resume_run_dir) if resume_run_dir else None,
            "bam_filtering": cfg["bam_filtering"],
            "cache": cfg["cache"],
            "output_compaction": cfg["output_compaction"],
            "config_path": str(config_path),
            "source_config_path": cfg.get("source_config_path"),
            "started_at": datetime.now().isoformat(timespec="seconds"),
        }
        with open(run_dir / "run_params.json", "w") as f:
            json.dump(run_params, f, indent=2)
    else:
        resume_event = {
            "event": "resume",
            "config_path": str(config_path),
            "resumed_at": datetime.now().isoformat(timespec="seconds"),
            "start_index": start_idx,
            "start_gene_id": int(cfg["gene_ids"][start_idx]),
            "start_gene_status": resume_plan["start_status"],
            "progress_snapshot": resume_plan["counts"],
        }
        _append_jsonl(run_dir / "resume_invocations.jsonl", resume_event)

    logger.info(f"Run directory: {run_dir}")
    if resume_mode:
        logger.info(f"Resuming existing run from gene index {start_idx} (gene_id={cfg['gene_ids'][start_idx]})")
        logger.info(f"Resume progress snapshot: {resume_plan['counts']}")
    elif resume_plan is not None:
        logger.info(
            "resume_run_dir provided but all genes are terminal (success={}, failed={}); starting new run",
            resume_plan["counts"]["success"],
            resume_plan["counts"]["failed"],
        )
    logger.info(f"Gene IDs: {cfg['gene_ids']}")
    logger.info(f"BLAST parameters: hitlist_size={cfg['hitlist_size']}, expect={cfg['blast_expect']}")
    logger.info(f"Read generation config: {cfg['read_generation']}")
    logger.info(f"Variant calling config: {cfg['variant_calling']}")
    logger.info(f"Ortholog selection config: {cfg['ortholog_selection']}")
    logger.info(f"BAM filtering config: {cfg['bam_filtering']}")
    logger.info(f"Cache config: {cfg['cache']}")
    logger.info(f"Output compaction config: {cfg['output_compaction']}")

    # Optional NCBI ortholog prefetch (batch) into cache.
    cache_cfg = cfg["cache"]
    batch_size = cache_cfg["ortholog_batch_size"]
    prefetch_source = None
    prefetch_enabled = (
        cache_cfg["enabled"]
        and cache_cfg["orthologs_ncbi"]
        and batch_size > 1
    )
    if prefetch_enabled:
        try:
            source = get_source("ncbi_datasets")
            if source.is_available() and hasattr(source, "prefetch_to_cache"):
                prefetch_source = source
            else:
                logger.warning("NCBI ortholog prefetch is disabled: ncbi_datasets source unavailable")
        except Exception as e:
            logger.warning(f"NCBI ortholog prefetch init failed: {e}")

    # Process each gene (ortholog source: NCBI Datasets -> BLAST fallback)
    success_count = 0
    failed_count = 0
    skipped_count = 0
    failed_genes_jsonl = run_dir / "failed_genes.jsonl"
    for idx in range(start_idx, len(cfg["gene_ids"])):
        gene_id = cfg["gene_ids"][idx]
        if prefetch_source is not None and ((idx - start_idx) % batch_size == 0):
            upcoming = cfg["gene_ids"][idx : idx + batch_size]
            try:
                logger.info(
                    "Prefetching NCBI orthologs for next {} genes",
                    len(upcoming),
                )
                produced = prefetch_source.prefetch_to_cache(
                    gene_ids=upcoming,
                    cache_dir=config.NCBI_ORTHOLOG_CACHE_DIR,
                    ortholog_scope=cfg["ortholog_selection"]["scope"],
                )
                logger.info(
                    "NCBI ortholog prefetch completed: {} cached gene(s)",
                    len(produced),
                )
            except Exception as e:
                logger.warning(f"NCBI ortholog prefetch failed: {e}")

        gene_dir = run_dir / f"gene_{gene_id}"
        status = _gene_status(run_dir, gene_id)
        if status == "success":
            logger.info(
                "Skipping gene {}: already has final artifact ({})",
                gene_id,
                GENE_SUCCESS_ARTIFACT,
            )
            skipped_count += 1
            continue
        if status == "failed":
            logger.info(
                "Skipping gene {}: already has failure record ({})",
                gene_id,
                GENE_FAILURE_ARTIFACT,
            )
            skipped_count += 1
            continue
        if status == "incomplete":
            logger.warning(
                "Gene {} has incomplete artifacts in {}; removing directory and rerunning from scratch",
                gene_id,
                gene_dir,
            )
            shutil.rmtree(gene_dir)

        started_at = datetime.now()
        try:
            run_gene(gene_id, gene_dir, cfg)
            success_count += 1
        except Exception as e:
            logger.exception(f"Failed to process gene {gene_id}")
            failed_count += 1
            finished_at = datetime.now()
            failure = {
                "gene_id": int(gene_id),
                "gene_dir": str(gene_dir),
                "status": "failed",
                "started_at": started_at.isoformat(timespec="seconds"),
                "finished_at": finished_at.isoformat(timespec="seconds"),
                "duration_seconds": round(
                (finished_at - started_at).total_seconds(),
                3,
                ),
                "error_type": type(e).__name__,
                "error_message": str(e) if str(e) else repr(e),
                "error_traceback": traceback.format_exc(),
            }
            gene_dir.mkdir(parents=True, exist_ok=True)
            _write_json(gene_dir / "failure.json", failure)
            _append_jsonl(failed_genes_jsonl, failure)

    logger.success(
        "Pipeline run complete (succeeded={}, failed={}, skipped={})",
        success_count,
        failed_count,
        skipped_count,
    )
    if failed_count > 0:
        logger.warning("Failure report: {}", failed_genes_jsonl)

    if cfg["output_compaction"]["enabled"]:
        logger.info("Output compaction enabled; writing compact statistics and deleting raw gene directories")
        compaction_result = compact_run_in_place(run_dir)
        logger.success(
            "Output compaction complete (success={}, failed={}, variants={}, deleted_gene_dirs={})",
            compaction_result.successful_genes,
            compaction_result.failed_genes,
            compaction_result.total_variants,
            compaction_result.deleted_gene_dirs,
        )


if __name__ == "__main__":
    main()
