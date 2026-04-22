from __future__ import annotations

import sys
from pathlib import Path
from typing import Iterable

import pandas as pd

ANALYSIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = ANALYSIS_DIR.parent
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from compare_runs_on_csv import build_outputs  # noqa: E402
from run_snapshot import (  # noqa: E402
    RunExportSummary,
    discover_exported_run_dirs,
    export_runs,
    load_exported_run_data,
)


def export_and_compare_runs(
    raw_runs_root: Path | str = REPO_ROOT / "runs",
    export_root: Path | str = REPO_ROOT / "analysis" / "data_runs",
    results_dir: Path | str = REPO_ROOT / "analysis" / "results" / "compact_compare",
    selected_runs: Iterable[str] | None = None,
    run_comparison: bool = True,
) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    """Export raw pipeline runs to compact snapshots and optionally compare them.

    This is the notebook-friendly wrapper for the analysis workflow. It does not
    delete raw run directories; use the returned table to decide what is safe to
    remove manually.
    """
    raw_runs_root = Path(raw_runs_root)
    export_root = Path(export_root)
    results_dir = Path(results_dir)
    selected_runs = list(selected_runs or [])

    summaries = export_runs(raw_runs_root, export_root, selected_runs=selected_runs)
    export_table = _build_export_table(summaries)

    outputs: dict[str, pd.DataFrame] = {}
    if run_comparison and summaries:
        run_data = {}
        exported_dirs = _select_exported_dirs(export_root, selected_runs)
        for export_dir in exported_dirs:
            run_label = export_dir.relative_to(export_root).as_posix()
            if run_label == ".":
                run_label = export_dir.name
            run_data[run_label] = load_exported_run_data(
                export_dir,
                run_label=run_label,
                include_variants=False,
            )
        if run_data:
            outputs = build_outputs(run_data, results_dir)

    return export_table, outputs


def _build_export_table(summaries: list[RunExportSummary]) -> pd.DataFrame:
    columns = [
        "run",
        "export_dir",
        "successful_genes",
        "failed_genes",
        "incomplete_genes",
        "missing_genes",
        "total_variants",
        "safe_to_delete_raw",
        "next_action",
    ]
    rows = []
    for summary in summaries:
        safe_to_delete_raw = summary.incomplete_genes == 0 and summary.missing_genes == 0
        rows.append(
            {
                "run": summary.relative_run_path,
                "export_dir": str(summary.output_dir),
                "successful_genes": summary.successful_genes,
                "failed_genes": summary.failed_genes,
                "incomplete_genes": summary.incomplete_genes,
                "missing_genes": summary.missing_genes,
                "total_variants": summary.total_variants,
                "safe_to_delete_raw": safe_to_delete_raw,
                "next_action": "review failures, then delete raw" if safe_to_delete_raw else "keep raw and resume/re-export",
            }
        )
    return pd.DataFrame(rows, columns=columns)


def _select_exported_dirs(export_root: Path, selected_runs: list[str]) -> list[Path]:
    run_dirs = discover_exported_run_dirs(export_root)
    if not selected_runs:
        return run_dirs

    selected = set(selected_runs)
    return [
        run_dir
        for run_dir in run_dirs
        if run_dir.name in selected or run_dir.relative_to(export_root).as_posix() in selected
    ]
