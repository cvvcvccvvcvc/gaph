# data_runs

This directory stores compact per-run analysis snapshots.

Recommended notebook workflow:

1. Open `analysis/export_and_compare_runs.ipynb`.
2. Set `RAW_RUNS_ROOT` to the folder that contains raw pipeline runs.
3. Press Run All.
4. Check the `safe_to_delete_raw` column.
5. Delete a raw run manually only when `safe_to_delete_raw=true` and the
   failed genes are expected/understood.

CLI equivalent, if needed:

1. Export raw runs into this directory:

```bash
python3 analysis/export_run_analysis_data.py runs/read_geometry_check --output-root analysis/data_runs/read_geometry_check
```

2. Run comparison on the exported snapshots:

```bash
python3 analysis/compare_runs_on_csv.py analysis/data_runs/read_geometry_check --output-dir analysis/results/read_geometry_check_csv
```

Per exported run, the snapshot contains:

- `run_manifest.json`: export metadata and row counts
- `run_params.json`: copied run parameters
- `genes.csv.gz`: one row per gene with final status and per-gene counts
- `variants.csv.gz`: compact per-variant table parsed from `gene_snps_annotated.vcf`
- `failure_events.csv.gz`: exported rows from `failed_genes.jsonl`
- `analysis_summary.json.gz`: precomputed run-level summary cache for fast notebook comparisons
- `ortholog_resolution.csv`: copied run-level ortholog audit, when present

After validating that the exported snapshot is sufficient for the analysis you need,
the original heavy run directory can be removed.

Before deleting a raw run, check the export command output or `run_manifest.json`.
`incomplete` and `missing` should be expected and understood; otherwise keep the raw
run and resume/re-export it first.
