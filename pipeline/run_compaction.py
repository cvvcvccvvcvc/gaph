from __future__ import annotations

import csv
import gzip
import json
import re
import shutil
import statistics
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator


EXPORT_SCHEMA_VERSION = 1

VARIANT_FIELDS = [
    "run_id",
    "gene_id",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_key",
    "filter",
    "clnsig",
    "clnrevstat",
    "clnsigconf",
    "af",
    "af_raw",
    "mc",
    "csq",
    "hgvsc",
    "has_clnsig",
    "has_clinvar_record",
    "has_gnomad",
    "has_both",
]

GENE_FIELDS = [
    "run_id",
    "gene_id",
    "status",
    "total_variants",
    "clinvar_intersected",
    "clinvar_record_count",
    "gnomad_intersected",
    "clinvar_gnomad_both",
    "mc_annotated_variants",
    "csq_annotated_variants",
    "error_type",
    "error_message",
]

FAILURE_EVENT_FIELDS = [
    "run_id",
    "gene_id",
    "status",
    "started_at",
    "finished_at",
    "duration_seconds",
    "error_type",
    "error_message",
    "gene_dir",
    "error_traceback",
]

_INFO_KEYS = {"CLNSIG", "CLNREVSTAT", "CLNSIGCONF", "AF", "MC", "CSQ", "HGVSC"}


@dataclass(frozen=True)
class RunExportSummary:
    run_label: str
    relative_run_path: str
    output_dir: Path
    successful_genes: int
    failed_genes: int
    incomplete_genes: int
    missing_genes: int
    total_variants: int


@dataclass(frozen=True)
class RunCompactionResult:
    run_dir: Path
    successful_genes: int
    failed_genes: int
    deleted_gene_dirs: int
    total_variants: int


COMPACT_RUN_FILES = [
    "run_manifest.json",
    "genes.csv.gz",
    "variants.csv.gz",
    "failure_events.csv.gz",
    "analysis_summary.json.gz",
]


def discover_raw_run_dirs(source_root: Path) -> list[Path]:
    source_root = Path(source_root)
    run_dirs = sorted({path.parent for path in source_root.rglob("run_params.json")})
    return [path for path in run_dirs if path.is_dir()]


def discover_exported_run_dirs(exports_root: Path) -> list[Path]:
    exports_root = Path(exports_root)
    run_dirs = sorted({path.parent for path in exports_root.rglob("run_manifest.json")})
    return [path for path in run_dirs if path.is_dir()]


def export_runs(
    source_root: Path,
    output_root: Path,
    selected_runs: Iterable[str] | None = None,
) -> list[RunExportSummary]:
    source_root = Path(source_root)
    output_root = Path(output_root)
    run_dirs = discover_raw_run_dirs(source_root)

    selected = {str(item) for item in selected_runs or []}
    if selected:
        run_dirs = [
            run_dir
            for run_dir in run_dirs
            if run_dir.name in selected or run_dir.relative_to(source_root).as_posix() in selected
        ]

    summaries = []
    for run_dir in run_dirs:
        relative_run_path = run_dir.relative_to(source_root)
        export_dir = output_root / relative_run_path
        summaries.append(export_run(run_dir, export_dir, relative_run_path.as_posix()))
    return summaries


def export_run(run_dir: Path, export_dir: Path, relative_run_path: str | None = None) -> RunExportSummary:
    run_dir = Path(run_dir)
    export_dir = Path(export_dir)
    export_dir.mkdir(parents=True, exist_ok=True)

    run_label = run_dir.name
    params = _load_run_params(run_dir)
    requested_gene_ids = [str(gene_id) for gene_id in params.get("gene_ids", []) if str(gene_id)]
    failure_events = _read_failure_events(run_dir, run_label)
    terminal_failures = _build_terminal_failures(run_dir, failure_events)
    failure_event_gene_ids = {str(event.get("gene_id", "")) for event in failure_events}
    for gene_id in sorted(set(terminal_failures) - failure_event_gene_ids, key=_gene_sort_key):
        failure_events.append(terminal_failures[gene_id])

    gene_rows: list[dict[str, Any]] = []
    seen_gene_ids: set[str] = set()
    summary_builder = AnalysisSummaryBuilder()
    variant_row_count = 0

    ortholog_resolution_path = run_dir / "ortholog_resolution.csv"
    exported_ortholog_resolution_path = export_dir / "ortholog_resolution.csv"
    if ortholog_resolution_path.exists():
        exported_ortholog_resolution_path.write_text(ortholog_resolution_path.read_text())
    else:
        exported_ortholog_resolution_path.unlink(missing_ok=True)

    _write_json(export_dir / "run_params.json", params)

    with gzip.open(export_dir / "variants.csv.gz", "wt", newline="") as handle:
        variant_writer = csv.DictWriter(handle, fieldnames=VARIANT_FIELDS)
        variant_writer.writeheader()

        for gene_dir in sorted(run_dir.glob("gene_*")):
            if not gene_dir.is_dir():
                continue

            gene_id = gene_dir.name.replace("gene_", "", 1)
            seen_gene_ids.add(gene_id)
            annotated_vcf = gene_dir / "gene_snps_annotated.vcf"

            if _is_nonempty_file(annotated_vcf):
                gene_stats = _new_gene_variant_stats()
                for row in iter_annotated_vcf_rows(annotated_vcf, run_label, gene_id):
                    variant_writer.writerow({key: row.get(key, "") for key in VARIANT_FIELDS})
                    summary_builder.add(row)
                    _update_gene_variant_stats(gene_stats, row)
                    variant_row_count += 1
                gene_rows.append(_build_gene_row_from_stats(run_label, gene_id, "success", gene_stats))
                continue

            failure_record = terminal_failures.get(gene_id, {})
            status = "failed" if failure_record else "incomplete"
            gene_rows.append(_build_gene_row(run_label, gene_id, status, [], failure_record))

    for gene_id in sorted(set(requested_gene_ids) | set(terminal_failures.keys())):
        if gene_id in seen_gene_ids:
            continue
        failure_record = terminal_failures.get(gene_id, {})
        status = "failed" if failure_record else "missing"
        gene_rows.append(_build_gene_row(run_label, gene_id, status, [], failure_record))

    _write_gzip_csv(export_dir / "genes.csv.gz", GENE_FIELDS, gene_rows)
    _write_gzip_csv(export_dir / "failure_events.csv.gz", FAILURE_EVENT_FIELDS, failure_events)
    _write_gzip_json(export_dir / "analysis_summary.json.gz", summary_builder.build())

    gene_counts = {
        "successful_genes": sum(1 for row in gene_rows if row["status"] == "success"),
        "failed_genes": sum(1 for row in gene_rows if row["status"] == "failed"),
        "incomplete_genes": sum(1 for row in gene_rows if row["status"] == "incomplete"),
        "missing_genes": sum(1 for row in gene_rows if row["status"] == "missing"),
    }

    manifest = {
        "schema_version": EXPORT_SCHEMA_VERSION,
        "exported_at": _utc_now_iso(),
        "run_label": run_label,
        "relative_run_path": relative_run_path or run_dir.name,
        "source_run_dir": str(run_dir.resolve()),
        "counts": {
            "requested_genes": len(requested_gene_ids),
            "gene_rows": len(gene_rows),
            "variant_rows": variant_row_count,
            "failure_events": len(failure_events),
            **gene_counts,
        },
        "files": {
            "run_params": "run_params.json",
            "genes": "genes.csv.gz",
            "variants": "variants.csv.gz",
            "failure_events": "failure_events.csv.gz",
            "analysis_summary": "analysis_summary.json.gz",
            "ortholog_resolution": "ortholog_resolution.csv" if ortholog_resolution_path.exists() else None,
        },
    }
    _write_json(export_dir / "run_manifest.json", manifest)

    return RunExportSummary(
        run_label=run_label,
        relative_run_path=manifest["relative_run_path"],
        output_dir=export_dir,
        successful_genes=gene_counts["successful_genes"],
        failed_genes=gene_counts["failed_genes"],
        incomplete_genes=gene_counts["incomplete_genes"],
        missing_genes=gene_counts["missing_genes"],
        total_variants=variant_row_count,
    )


def compact_run_in_place(run_dir: Path) -> RunCompactionResult:
    """Write compact run statistics into run_dir and remove raw per-gene outputs.

    The destructive part is intentionally last:
    1. export compact files into a temporary directory;
    2. validate the temporary snapshot;
    3. publish compact files into the run directory;
    4. validate the published snapshot;
    5. delete raw gene_* directories.
    """
    run_dir = Path(run_dir)
    if not run_dir.exists() or not run_dir.is_dir():
        raise FileNotFoundError(f"Run directory does not exist: {run_dir}")

    with tempfile.TemporaryDirectory(prefix=".compaction_", dir=run_dir) as tmp_name:
        tmp_dir = Path(tmp_name)
        summary = export_run(run_dir, tmp_dir, run_dir.name)
        validate_compacted_run(tmp_dir, allow_incomplete=False)

        for filename in COMPACT_RUN_FILES:
            source = tmp_dir / filename
            if not source.exists():
                raise FileNotFoundError(f"Compaction did not produce required file: {source}")
            source.replace(run_dir / filename)

    published = validate_compacted_run(run_dir, allow_incomplete=False)
    _mark_compaction(
        run_dir,
        raw_gene_dirs_deleted=False,
        deleted_gene_dirs=0,
        validation=published,
    )
    validate_compacted_run(run_dir, allow_incomplete=False)

    deleted_gene_dirs = _delete_raw_gene_dirs(run_dir)
    _mark_compaction(
        run_dir,
        raw_gene_dirs_deleted=True,
        deleted_gene_dirs=deleted_gene_dirs,
        validation=published,
    )
    validate_compacted_run(run_dir, allow_incomplete=False)
    _delete_compaction_temp_dirs(run_dir)

    return RunCompactionResult(
        run_dir=run_dir,
        successful_genes=summary.successful_genes,
        failed_genes=summary.failed_genes,
        deleted_gene_dirs=deleted_gene_dirs,
        total_variants=summary.total_variants,
    )


def validate_compacted_run(export_dir: Path, *, allow_incomplete: bool = False) -> dict[str, Any]:
    """Validate that a compact snapshot is readable and internally consistent."""
    export_dir = Path(export_dir)
    manifest_path = export_dir / "run_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(f"Missing compact run manifest: {manifest_path}")

    for filename in COMPACT_RUN_FILES[1:]:
        path = export_dir / filename
        if not _is_nonempty_file(path):
            raise FileNotFoundError(f"Missing or empty compact run file: {path}")

    manifest = _load_json(manifest_path)
    counts = manifest.get("counts", {})
    if not isinstance(counts, dict):
        raise ValueError(f"Invalid compact manifest counts: {manifest_path}")

    genes_rows = _read_gzip_csv(export_dir / "genes.csv.gz")
    variant_row_count = _count_gzip_csv_rows(export_dir / "variants.csv.gz")
    failure_rows = _read_gzip_csv(export_dir / "failure_events.csv.gz")
    _read_gzip_json(export_dir / "analysis_summary.json.gz")

    expected_counts = {
        "gene_rows": len(genes_rows),
        "variant_rows": variant_row_count,
        "failure_events": len(failure_rows),
        "successful_genes": sum(1 for row in genes_rows if row.get("status") == "success"),
        "failed_genes": sum(1 for row in genes_rows if row.get("status") == "failed"),
        "incomplete_genes": sum(1 for row in genes_rows if row.get("status") == "incomplete"),
        "missing_genes": sum(1 for row in genes_rows if row.get("status") == "missing"),
    }
    for key, actual in expected_counts.items():
        expected = _as_int(counts.get(key, 0))
        if expected != actual:
            raise ValueError(
                f"Compact run manifest count mismatch for {key}: "
                f"manifest={expected}, actual={actual}"
            )

    if not allow_incomplete and (
        expected_counts["incomplete_genes"] > 0 or expected_counts["missing_genes"] > 0
    ):
        raise ValueError(
            "Refusing to delete raw outputs because compact snapshot contains "
            f"incomplete={expected_counts['incomplete_genes']} and "
            f"missing={expected_counts['missing_genes']} gene(s)"
        )

    requested = {
        str(gene_id)
        for gene_id in _load_run_params(export_dir).get("gene_ids", [])
        if str(gene_id)
    }
    observed = {str(row.get("gene_id", "")) for row in genes_rows if row.get("gene_id", "")}
    if requested and not requested.issubset(observed):
        missing = sorted(requested - observed, key=_gene_sort_key)
        raise ValueError(f"Compact snapshot is missing requested gene row(s): {missing}")

    return {
        "counts": expected_counts,
        "manifest": manifest,
    }


def load_compacted_gene_status(run_dir: Path, gene_id: int | str) -> str | None:
    """Return terminal status from compacted genes.csv.gz when raw gene files are gone."""
    run_dir = Path(run_dir)
    if not (run_dir / "run_manifest.json").exists() or not (run_dir / "genes.csv.gz").exists():
        return None

    wanted = str(gene_id)
    for row in _read_gzip_csv(run_dir / "genes.csv.gz"):
        if str(row.get("gene_id", "")) == wanted:
            status = str(row.get("status", "")).strip()
            return status or None
    return None


def run_has_compacted_snapshot(run_dir: Path) -> bool:
    run_dir = Path(run_dir)
    return (run_dir / "run_manifest.json").exists() and (run_dir / "genes.csv.gz").exists()


def _delete_raw_gene_dirs(run_dir: Path) -> int:
    deleted = 0
    for path in sorted(run_dir.glob("gene_*")):
        if not path.is_dir():
            continue
        shutil.rmtree(path)
        deleted += 1
    return deleted


def _delete_compaction_temp_dirs(run_dir: Path) -> int:
    deleted = 0
    for path in sorted(run_dir.glob(".compaction_*")):
        if not path.is_dir():
            continue
        shutil.rmtree(path)
        deleted += 1
    return deleted


def _mark_compaction(
    run_dir: Path,
    *,
    raw_gene_dirs_deleted: bool,
    deleted_gene_dirs: int,
    validation: dict[str, Any],
) -> None:
    manifest_path = run_dir / "run_manifest.json"
    manifest = _load_json(manifest_path)
    manifest["compaction"] = {
        "schema_version": 1,
        "compacted_at": _utc_now_iso(),
        "raw_gene_dirs_deleted": raw_gene_dirs_deleted,
        "deleted_gene_dirs": deleted_gene_dirs,
        "successful_genes": validation["counts"]["successful_genes"],
        "failed_genes": validation["counts"]["failed_genes"],
        "total_variants": validation["counts"]["variant_rows"],
    }
    _write_json(manifest_path, manifest)


def iter_annotated_vcf_rows(
    vcf_path: Path,
    run_label: str,
    gene_id: str,
) -> Iterator[dict[str, Any]]:
    with Path(vcf_path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue

            chrom, pos, _variant_id, ref, alt, _qual, filter_value, info = parts[:8]
            info_values = _extract_info_values(info)
            clnsig = info_values.get("CLNSIG", "")
            clnrevstat = info_values.get("CLNREVSTAT", "")
            clnsigconf = info_values.get("CLNSIGCONF", "")
            af_raw = info_values.get("AF", "")
            af = _parse_float(af_raw)
            mc = info_values.get("MC", "")
            csq = info_values.get("CSQ", "")
            hgvsc = info_values.get("HGVSC", "")

            has_clnsig = bool(clnsig)
            has_clinvar_record = bool(clnsig or clnrevstat or clnsigconf)
            has_gnomad = af is not None

            yield {
                "run_id": run_label,
                "gene_id": str(gene_id),
                "chrom": chrom,
                "pos": int(pos),
                "ref": ref,
                "alt": alt,
                "variant_key": f"{chrom}:{pos}:{ref}>{alt}",
                "filter": filter_value,
                "clnsig": clnsig,
                "clnrevstat": clnrevstat,
                "clnsigconf": clnsigconf,
                "af": af if af is not None else "",
                "af_raw": af_raw,
                "mc": mc,
                "csq": csq,
                "hgvsc": hgvsc,
                "has_clnsig": int(has_clnsig),
                "has_clinvar_record": int(has_clinvar_record),
                "has_gnomad": int(has_gnomad),
                "has_both": int(has_clnsig and has_gnomad),
            }


def parse_annotated_vcf(vcf_path: Path, run_label: str, gene_id: str) -> list[dict[str, Any]]:
    return list(iter_annotated_vcf_rows(vcf_path, run_label, gene_id))


def load_exported_run_data(
    export_dir: Path,
    genes_filter: set[str] | None = None,
    run_label: str | None = None,
    include_variants: bool = True,
) -> dict[str, Any]:
    export_dir = Path(export_dir)
    params = _load_json(export_dir / "run_params.json")
    manifest = _load_json(export_dir / "run_manifest.json")
    genes_rows = _read_gzip_csv(export_dir / "genes.csv.gz")
    failure_events = _read_gzip_csv(export_dir / "failure_events.csv.gz")
    analysis_summary = None

    if genes_filter:
        genes_rows = [row for row in genes_rows if str(row.get("gene_id", "")) in genes_filter]
        failure_events = [row for row in failure_events if str(row.get("gene_id", "")) in genes_filter]
    else:
        analysis_summary = load_exported_analysis_summary(export_dir)

    data = {
        "run_id": run_label or manifest.get("relative_run_path") or manifest.get("run_label") or export_dir.name,
        "params": params,
        "manifest": manifest,
        "gene_rows": genes_rows,
        "variant_rows": [],
        "failure_events": failure_events,
        "gene_results": _build_gene_results_from_gene_rows(genes_rows),
        "all_clnsig": [],
        "all_clinvar_records": [],
        "all_gnomad_af": [],
        "analysis_summary": analysis_summary,
    }
    if include_variants:
        return hydrate_exported_run_data(data, export_dir, genes_filter=genes_filter)
    return data


def hydrate_exported_run_data(
    data: dict[str, Any],
    export_dir: Path,
    genes_filter: set[str] | None = None,
) -> dict[str, Any]:
    export_dir = Path(export_dir)
    variant_rows = _read_gzip_csv(export_dir / "variants.csv.gz")

    if genes_filter:
        variant_rows = [row for row in variant_rows if str(row.get("gene_id", "")) in genes_filter]

    variants_by_gene: dict[str, list[dict[str, str]]] = {}
    for row in variant_rows:
        gene_id = str(row.get("gene_id", ""))
        variants_by_gene.setdefault(gene_id, []).append(row)

    all_clnsig = [row["clnsig"] for row in variant_rows if row.get("clnsig")]
    all_clinvar_records = [
        {
            "clnsig": row.get("clnsig", ""),
            "clnrevstat": row.get("clnrevstat", ""),
            "clnsigconf": row.get("clnsigconf", ""),
        }
        for row in variant_rows
        if _as_int(row.get("has_clinvar_record", 0))
    ]
    all_gnomad_af = [
        float(row["af"])
        for row in variant_rows
        if str(row.get("af", "")).strip() != ""
    ]

    hydrated_gene_results = []
    for gene_result in data.get("gene_results", []):
        gene_id = str(gene_result["gene_id"])
        gene_variants = variants_by_gene.get(gene_id, [])
        hydrated_gene_results.append(
            {
                **gene_result,
                "clnsig_values": [variant["clnsig"] for variant in gene_variants if variant.get("clnsig")],
                "clinvar_records": [
                    {
                        "clnsig": variant.get("clnsig", ""),
                        "clnrevstat": variant.get("clnrevstat", ""),
                        "clnsigconf": variant.get("clnsigconf", ""),
                    }
                    for variant in gene_variants
                    if _as_int(variant.get("has_clinvar_record", 0))
                ],
                "gnomad_af_values": [
                    float(variant["af"])
                    for variant in gene_variants
                    if str(variant.get("af", "")).strip() != ""
                ],
            }
        )

    return {
        **data,
        "variant_rows": variant_rows,
        "gene_results": hydrated_gene_results,
        "all_clnsig": all_clnsig,
        "all_clinvar_records": all_clinvar_records,
        "all_gnomad_af": all_gnomad_af,
        "analysis_summary": data.get("analysis_summary"),
    }


def load_exported_analysis_summary(export_dir: Path) -> dict[str, Any] | None:
    export_dir = Path(export_dir)
    summary_path = export_dir / "analysis_summary.json.gz"
    if summary_path.exists():
        return _read_gzip_json(summary_path)

    variants_path = export_dir / "variants.csv.gz"
    if not variants_path.exists():
        return None

    summary_builder = AnalysisSummaryBuilder()
    for row in _iter_gzip_csv(variants_path):
        summary_builder.add(row)
    summary = summary_builder.build()
    _write_gzip_json(summary_path, summary)
    return summary


def build_analysis_summary_from_variant_rows(variant_rows: list[dict[str, Any]]) -> dict[str, Any]:
    summary_builder = AnalysisSummaryBuilder()
    for row in variant_rows:
        summary_builder.add(row)
    return summary_builder.build()


class AnalysisSummaryBuilder:
    """Incrementally build run-level analysis without retaining every variant row."""

    def __init__(self) -> None:
        self.clnsig_counts: dict[str, int] = {}
        self.grouped_clnsig_counts: dict[str, int] = {}
        self.benign_star_counts: dict[str, int] = {}
        self.conflicting_submitter_counts: dict[str, int] = {}
        self.conflicting_gene_counts: dict[str, dict[str, int]] = {}
        self.mc_term_counts: dict[str, int] = {}
        self.csq_term_counts: dict[str, int] = {}
        self.af_values: list[float] = []
        self.both_annotated_count = 0
        self.discordant_count = 0
        self.gene_both_counts: dict[str, int] = {}
        self.gene_discordant_counts: dict[str, int] = {}
        self.impact_run_counts: dict[tuple[str, str], int] = {}
        self.impact_run_source_totals: dict[str, int] = {}
        self.impact_gene_counts: dict[tuple[str, str, str], int] = {}
        self.impact_gene_source_totals: dict[tuple[str, str], int] = {}

    def add(self, row: dict[str, Any]) -> None:
        gene_id = str(row.get("gene_id", ""))
        clnsig = str(row.get("clnsig", "")).strip()
        clnrevstat = str(row.get("clnrevstat", "")).strip()
        clnsigconf = str(row.get("clnsigconf", "")).strip()

        if clnsig:
            for label in split_clnsig_labels(clnsig):
                self.clnsig_counts[label] = self.clnsig_counts.get(label, 0) + 1
            grouped = group_clnsig(clnsig)
            self.grouped_clnsig_counts[grouped] = self.grouped_clnsig_counts.get(grouped, 0) + 1

            labels = split_clnsig_labels(clnsig)
            if any(is_benign_label(label) for label in labels):
                stars = clinvar_stars_from_review_status(clnrevstat)
                bucket = f"benign_{stars}" if stars is not None else "benign_unknown"
                self.benign_star_counts[bucket] = self.benign_star_counts.get(bucket, 0) + 1

            if "Conflicting_classifications_of_pathogenicity" in labels:
                parsed = parse_clnsigconf_submitter_counts(clnsigconf)
                gene_bucket = self.conflicting_gene_counts.setdefault(gene_id, {})
                for label, count in parsed.items():
                    self.conflicting_submitter_counts[label] = self.conflicting_submitter_counts.get(label, 0) + count
                    gene_bucket[label] = gene_bucket.get(label, 0) + count

        af_raw = str(row.get("af", "")).strip()
        if af_raw:
            self.af_values.append(float(af_raw))

        mc_terms = parse_mc_terms(str(row.get("mc", "")))
        csq_terms = parse_csq_terms(str(row.get("csq", "")))
        mc_class = classify_terms(mc_terms)
        csq_class = classify_terms(csq_terms)

        for term in mc_terms:
            self.mc_term_counts[term] = self.mc_term_counts.get(term, 0) + 1
        for term in csq_terms:
            self.csq_term_counts[term] = self.csq_term_counts.get(term, 0) + 1

        if mc_class is not None:
            self._add_impact(gene_id, "ClinVar_MC", mc_class)
        if csq_class is not None:
            self._add_impact(gene_id, "gnomAD_CSQ", csq_class)

        if mc_class is not None and csq_class is not None:
            self.both_annotated_count += 1
            self.gene_both_counts[gene_id] = self.gene_both_counts.get(gene_id, 0) + 1
            if mc_class != csq_class:
                self.discordant_count += 1
                self.gene_discordant_counts[gene_id] = self.gene_discordant_counts.get(gene_id, 0) + 1

    def _add_impact(self, gene_id: str, source: str, impact_class: str) -> None:
        self.impact_run_counts[(source, impact_class)] = (
            self.impact_run_counts.get((source, impact_class), 0) + 1
        )
        self.impact_run_source_totals[source] = self.impact_run_source_totals.get(source, 0) + 1
        self.impact_gene_counts[(gene_id, source, impact_class)] = (
            self.impact_gene_counts.get((gene_id, source, impact_class), 0) + 1
        )
        self.impact_gene_source_totals[(gene_id, source)] = (
            self.impact_gene_source_totals.get((gene_id, source), 0) + 1
        )

    def build(self) -> dict[str, Any]:
        run_impact_rows, run_gene_impact_rows = self._summarize_impact_counts()
        run_mismatch_summary = (
            {
                "discordant_count": self.discordant_count,
                "both_annotated_count": self.both_annotated_count,
                "discordant_share_%": round(
                    self.discordant_count / self.both_annotated_count * 100,
                    2,
                )
                if self.both_annotated_count
                else 0.0,
            }
            if self.both_annotated_count
            else None
        )
        run_gene_mismatch_rows = [
            {
                "gene_id": gene_id,
                "discordant_count": self.gene_discordant_counts.get(gene_id, 0),
                "both_annotated_count": self.gene_both_counts.get(gene_id, 0),
                "discordant_share_%": round(
                    self.gene_discordant_counts.get(gene_id, 0)
                    / self.gene_both_counts.get(gene_id, 0)
                    * 100,
                    2,
                )
                if self.gene_both_counts.get(gene_id, 0)
                else 0.0,
            }
            for gene_id in sorted(self.gene_discordant_counts.keys(), key=_gene_sort_key)
        ]

        conflicting_gene_rows = []
        all_conf_labels = sorted(self.conflicting_submitter_counts.keys(), key=str.lower)
        for gene_id in sorted(self.conflicting_gene_counts.keys(), key=_gene_sort_key):
            row = {"gene_id": gene_id}
            per_gene = self.conflicting_gene_counts[gene_id]
            for label in all_conf_labels:
                row[label] = int(per_gene.get(label, 0))
            row["total_submitters"] = int(sum(per_gene.values()))
            conflicting_gene_rows.append(row)

        return {
            "schema_version": EXPORT_SCHEMA_VERSION,
            "clnsig_counts": _counts_to_rows("CLNSIG", self.clnsig_counts),
            "grouped_clnsig_counts": _counts_to_rows("Group", self.grouped_clnsig_counts),
            "benign_star_counts": _counts_to_rows("benign_star_bucket", self.benign_star_counts),
            "conflicting_clnsigconf_counts": _counts_to_rows(
                "CLNSIGCONF",
                self.conflicting_submitter_counts,
            ),
            "conflicting_clnsigconf_by_gene": conflicting_gene_rows,
            "af_stats": build_af_summary(self.af_values),
            "run_impact_rows": run_impact_rows,
            "run_gene_impact_rows": run_gene_impact_rows,
            "run_mismatch_summary": run_mismatch_summary,
            "run_gene_mismatch_rows": run_gene_mismatch_rows,
            "mc_terms": _counts_to_rows("MC_term", self.mc_term_counts),
            "csq_terms": _counts_to_rows("CSQ_term", self.csq_term_counts),
        }

    def _summarize_impact_counts(self) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
        run_rows = []
        for (source, impact_class), count in sorted(self.impact_run_counts.items()):
            source_total = self.impact_run_source_totals[source]
            run_rows.append(
                {
                    "source": source,
                    "impact_class": impact_class,
                    "count": count,
                    "source_total": source_total,
                    "share_%": round(count / source_total * 100, 2) if source_total else 0.0,
                }
            )

        gene_rows = []
        for (gene_id, source, impact_class), count in sorted(
            self.impact_gene_counts.items(),
            key=lambda item: (item[0][1], _gene_sort_key(item[0][0]), item[0][2]),
        ):
            source_total = self.impact_gene_source_totals[(gene_id, source)]
            gene_rows.append(
                {
                    "gene_id": gene_id,
                    "source": source,
                    "impact_class": impact_class,
                    "count": count,
                    "source_total": source_total,
                    "share_%": round(count / source_total * 100, 2) if source_total else 0.0,
                }
            )

        return run_rows, gene_rows


def _build_terminal_failures(run_dir: Path, failure_events: list[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    terminal: dict[str, dict[str, Any]] = {}

    for event in failure_events:
        gene_id = str(event.get("gene_id", "")).strip()
        if gene_id:
            terminal[gene_id] = event

    for gene_dir in sorted(run_dir.glob("gene_*")):
        failure_json = gene_dir / "failure.json"
        if not _is_nonempty_file(failure_json):
            continue
        gene_id = gene_dir.name.replace("gene_", "", 1)
        try:
            record = _load_json(failure_json)
        except json.JSONDecodeError as exc:
            terminal[gene_id] = {
                "run_id": run_dir.name,
                "gene_id": gene_id,
                "status": "failed",
                "started_at": "",
                "finished_at": "",
                "duration_seconds": "",
                "error_type": "InvalidFailureRecord",
                "error_message": f"Could not parse failure.json: {exc}",
                "gene_dir": str(gene_dir),
                "error_traceback": "",
            }
            continue

        gene_id = str(record.get("gene_id") or gene_id)
        terminal[gene_id] = {
            "run_id": run_dir.name,
            "gene_id": gene_id,
            "status": record.get("status", "failed"),
            "started_at": record.get("started_at", ""),
            "finished_at": record.get("finished_at", ""),
            "duration_seconds": record.get("duration_seconds", ""),
            "error_type": record.get("error_type", ""),
            "error_message": record.get("error_message", ""),
            "gene_dir": record.get("gene_dir", str(gene_dir)),
            "error_traceback": record.get("error_traceback", ""),
        }

    return terminal


def _read_failure_events(run_dir: Path, run_label: str) -> list[dict[str, Any]]:
    path = run_dir / "failed_genes.jsonl"
    if not path.exists():
        return []

    rows = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            payload = json.loads(line)
            rows.append(
                {
                    "run_id": run_label,
                    "gene_id": str(payload.get("gene_id", "")),
                    "status": payload.get("status", "failed"),
                    "started_at": payload.get("started_at", ""),
                    "finished_at": payload.get("finished_at", ""),
                    "duration_seconds": payload.get("duration_seconds", ""),
                    "error_type": payload.get("error_type", ""),
                    "error_message": payload.get("error_message", ""),
                    "gene_dir": payload.get("gene_dir", ""),
                    "error_traceback": payload.get("error_traceback", ""),
                }
            )

    return rows


def _build_gene_row(
    run_label: str,
    gene_id: str,
    status: str,
    variant_rows: list[dict[str, Any]],
    failure_record: dict[str, Any] | None = None,
) -> dict[str, Any]:
    stats = _new_gene_variant_stats()
    for row in variant_rows:
        _update_gene_variant_stats(stats, row)
    return _build_gene_row_from_stats(run_label, gene_id, status, stats, failure_record)


def _new_gene_variant_stats() -> dict[str, int]:
    return {
        "total_variants": 0,
        "clinvar_intersected": 0,
        "clinvar_record_count": 0,
        "gnomad_intersected": 0,
        "clinvar_gnomad_both": 0,
        "mc_annotated_variants": 0,
        "csq_annotated_variants": 0,
    }


def _update_gene_variant_stats(stats: dict[str, int], row: dict[str, Any]) -> None:
    stats["total_variants"] += 1
    stats["clinvar_intersected"] += _as_int(row.get("has_clnsig", 0))
    stats["clinvar_record_count"] += _as_int(row.get("has_clinvar_record", 0))
    stats["gnomad_intersected"] += _as_int(row.get("has_gnomad", 0))
    stats["clinvar_gnomad_both"] += _as_int(row.get("has_both", 0))
    stats["mc_annotated_variants"] += int(bool(row.get("mc")))
    stats["csq_annotated_variants"] += int(bool(row.get("csq")))


def _build_gene_row_from_stats(
    run_label: str,
    gene_id: str,
    status: str,
    stats: dict[str, int],
    failure_record: dict[str, Any] | None = None,
) -> dict[str, Any]:
    failure_record = failure_record or {}
    return {
        "run_id": run_label,
        "gene_id": str(gene_id),
        "status": status,
        "total_variants": stats["total_variants"],
        "clinvar_intersected": stats["clinvar_intersected"],
        "clinvar_record_count": stats["clinvar_record_count"],
        "gnomad_intersected": stats["gnomad_intersected"],
        "clinvar_gnomad_both": stats["clinvar_gnomad_both"],
        "mc_annotated_variants": stats["mc_annotated_variants"],
        "csq_annotated_variants": stats["csq_annotated_variants"],
        "error_type": failure_record.get("error_type", ""),
        "error_message": failure_record.get("error_message", ""),
    }


def _build_gene_results_from_gene_rows(genes_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    successful_gene_rows = [row for row in genes_rows if row.get("status") == "success"]
    return [
        {
            "gene_id": str(row["gene_id"]),
            "total_variants": _as_int(row.get("total_variants", 0)),
            "clinvar_intersected": _as_int(row.get("clinvar_intersected", 0)),
            "gnomad_intersected": _as_int(row.get("gnomad_intersected", 0)),
            "clinvar_gnomad_both": _as_int(row.get("clinvar_gnomad_both", 0)),
            "clnsig_values": [],
            "clinvar_records": [],
            "gnomad_af_values": [],
        }
        for row in successful_gene_rows
    ]


def split_clnsig_labels(clnsig_raw: str) -> list[str]:
    if not clnsig_raw:
        return []
    return [item.strip() for item in re.split(r"[|,]", str(clnsig_raw)) if item.strip()]


def clinvar_stars_from_review_status(review_status: str) -> int | None:
    if not review_status:
        return None

    star_values = []
    for token in str(review_status).split("|"):
        token = token.strip().lower()
        if not token:
            continue
        if "practice_guideline" in token:
            star_values.append(4)
        elif "reviewed_by_expert_panel" in token:
            star_values.append(3)
        elif "criteria_provided,_multiple_submitters,_no_conflicts" in token:
            star_values.append(2)
        elif "criteria_provided,_single_submitter" in token or "criteria_provided,_conflicting_classifications" in token:
            star_values.append(1)
        elif (
            "no_assertion_criteria_provided" in token
            or "no_assertion_provided" in token
            or "no_interpretation_for_the_single_variant" in token
        ):
            star_values.append(0)
    return max(star_values) if star_values else None


def is_benign_label(label: str) -> bool:
    normalized = str(label).strip().lower()
    return ("benign" in normalized) and ("pathogenic" not in normalized) and ("conflicting" not in normalized)


def parse_clnsigconf_submitter_counts(clnsigconf_raw: str) -> dict[str, int]:
    counts: dict[str, int] = {}
    if not clnsigconf_raw:
        return counts
    for token in str(clnsigconf_raw).split("|"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^(.*)\((\d+)\)$", token)
        if match:
            label = match.group(1).strip()
            submitter_count = int(match.group(2))
        else:
            label = token
            submitter_count = 1
        if label:
            counts[label] = counts.get(label, 0) + submitter_count
    return counts


def group_clnsig(clnsig: str) -> str:
    clnsig_lower = str(clnsig).lower()
    if "pathogenic" in clnsig_lower and "benign" not in clnsig_lower and "conflicting" not in clnsig_lower:
        return "Pathogenic"
    if "benign" in clnsig_lower and "pathogenic" not in clnsig_lower and "conflicting" not in clnsig_lower:
        return "Benign"
    if "uncertain" in clnsig_lower:
        return "VUS"
    if "conflicting" in clnsig_lower or "not_provided" in clnsig_lower:
        return "Conflicting"
    return "Other"


IMPACT_RANK = {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3, "OTHER": -1}

TERM_TO_IMPACT = {
    "frameshift_variant": "HIGH",
    "stop_gained": "HIGH",
    "stop_lost": "HIGH",
    "start_lost": "HIGH",
    "splice_acceptor_variant": "HIGH",
    "splice_donor_variant": "HIGH",
    "missense_variant": "MODERATE",
    "inframe_insertion": "MODERATE",
    "inframe_deletion": "MODERATE",
    "protein_altering_variant": "MODERATE",
    "splice_region_variant": "MODERATE",
    "synonymous_variant": "LOW",
    "intron_variant": "MODIFIER",
    "5_prime_utr_variant": "MODIFIER",
    "3_prime_utr_variant": "MODIFIER",
    "upstream_gene_variant": "MODIFIER",
    "downstream_gene_variant": "MODIFIER",
    "intergenic_variant": "MODIFIER",
    "regulatory_region_variant": "MODIFIER",
    "non_coding_transcript_variant": "MODIFIER",
    "coding_sequence_variant": "MODIFIER",
    "mature_mirna_variant": "MODIFIER",
}


def parse_mc_terms(raw: str) -> list[str]:
    out = []
    if not raw:
        return out
    for item in str(raw).split(","):
        item = item.strip()
        if not item:
            continue
        if "|" in item:
            term = item.split("|", 1)[1].strip()
        else:
            term = item
        if term:
            out.append(term)
    return out


def parse_csq_terms(raw: str) -> list[str]:
    out = []
    if not raw:
        return out
    for item in re.split(r"[|,]", str(raw)):
        item = item.strip()
        if item:
            out.append(item)
    return out


def classify_terms(terms: list[str]) -> str | None:
    if not terms:
        return None
    best = "OTHER"
    best_rank = IMPACT_RANK[best]
    for term in terms:
        impact = TERM_TO_IMPACT.get(term, "OTHER")
        rank = IMPACT_RANK[impact]
        if rank > best_rank:
            best = impact
            best_rank = rank
    return best


def _summarize_impact_rows(impact_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    run_counts: dict[tuple[str, str], int] = {}
    run_source_totals: dict[str, int] = {}
    gene_counts: dict[tuple[str, str, str], int] = {}
    gene_source_totals: dict[tuple[str, str], int] = {}

    for row in impact_rows:
        source = str(row["source"])
        gene_id = str(row["gene_id"])
        impact_class = str(row["impact_class"])

        run_counts[(source, impact_class)] = run_counts.get((source, impact_class), 0) + 1
        run_source_totals[source] = run_source_totals.get(source, 0) + 1

        gene_counts[(gene_id, source, impact_class)] = gene_counts.get((gene_id, source, impact_class), 0) + 1
        gene_source_totals[(gene_id, source)] = gene_source_totals.get((gene_id, source), 0) + 1

    run_rows = []
    for (source, impact_class), count in sorted(run_counts.items()):
        source_total = run_source_totals[source]
        run_rows.append(
            {
                "source": source,
                "impact_class": impact_class,
                "count": count,
                "source_total": source_total,
                "share_%": round(count / source_total * 100, 2) if source_total else 0.0,
            }
        )

    gene_rows = []
    for (gene_id, source, impact_class), count in sorted(gene_counts.items(), key=lambda item: (item[0][1], _gene_sort_key(item[0][0]), item[0][2])):
        source_total = gene_source_totals[(gene_id, source)]
        gene_rows.append(
            {
                "gene_id": gene_id,
                "source": source,
                "impact_class": impact_class,
                "count": count,
                "source_total": source_total,
                "share_%": round(count / source_total * 100, 2) if source_total else 0.0,
            }
        )

    return run_rows, gene_rows


def build_af_summary(af_values: list[float]) -> dict[str, Any]:
    if not af_values:
        return {
            "count": 0,
            "min": None,
            "max": None,
            "mean": None,
            "median": None,
            "rare": 0,
            "low_freq": 0,
            "common": 0,
        }

    rare = sum(1 for value in af_values if value < 0.01)
    low_freq = sum(1 for value in af_values if 0.01 <= value < 0.05)
    common = sum(1 for value in af_values if value >= 0.05)

    return {
        "count": len(af_values),
        "min": min(af_values),
        "max": max(af_values),
        "mean": sum(af_values) / len(af_values),
        "median": statistics.median(af_values),
        "rare": rare,
        "low_freq": low_freq,
        "common": common,
    }


def _counts_to_rows(key_name: str, counts: dict[str, int]) -> list[dict[str, Any]]:
    return [{key_name: key, "count": int(counts[key])} for key in sorted(counts.keys(), key=str.lower)]


def _gene_sort_key(gene_id: str) -> tuple[int, str]:
    text = str(gene_id)
    return (0, f"{int(text):020d}") if text.isdigit() else (1, text)


def _extract_info_values(info: str) -> dict[str, str]:
    values: dict[str, str] = {}
    for token in info.split(";"):
        if "=" not in token:
            continue
        key, value = token.split("=", 1)
        if key in _INFO_KEYS:
            values[key] = value
    return values


def _parse_float(value: str) -> float | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        return None


def _is_nonempty_file(path: Path) -> bool:
    return path.exists() and path.is_file() and path.stat().st_size > 0


def _as_int(value: Any) -> int:
    if value is None:
        return 0
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, int):
        return value
    text = str(value).strip()
    if not text:
        return 0
    return int(float(text))


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")


def _write_gzip_json(path: Path, payload: dict[str, Any]) -> None:
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle, ensure_ascii=False)


def _write_gzip_csv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def _read_gzip_json(path: Path) -> dict[str, Any]:
    with gzip.open(path, "rt") as handle:
        return json.load(handle)


def _read_gzip_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle))


def _iter_gzip_csv(path: Path) -> Iterator[dict[str, str]]:
    if not path.exists():
        return
    with gzip.open(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle)


def _count_gzip_csv_rows(path: Path) -> int:
    if not path.exists():
        return 0
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _row in reader)


def _load_run_params(run_dir: Path) -> dict[str, Any]:
    path = run_dir / "run_params.json"
    return _load_json(path) if path.exists() else {}


def _load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text())


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()
