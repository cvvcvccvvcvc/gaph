from __future__ import annotations

import argparse
from pathlib import Path

from run_snapshot import export_runs


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Export compact per-run analysis snapshots from raw pipeline runs."
    )
    parser.add_argument(
        "source_root",
        type=Path,
        help="Root directory that contains raw run folders.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("analysis/data_runs"),
        help="Where compact run snapshots will be written.",
    )
    parser.add_argument(
        "--run",
        action="append",
        default=[],
        help="Run folder name or relative path under source_root. Can be passed multiple times.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    summaries = export_runs(args.source_root, args.output_root, selected_runs=args.run)
    if not summaries:
        print(f"No runs found under {args.source_root}")
        return 1

    print(f"Exported {len(summaries)} run snapshots to {args.output_root}")
    for summary in summaries:
        print(
            f"- {summary.relative_run_path}: "
            f"success={summary.successful_genes}, "
            f"failed={summary.failed_genes}, "
            f"incomplete={summary.incomplete_genes}, "
            f"missing={summary.missing_genes}, "
            f"variants={summary.total_variants}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
