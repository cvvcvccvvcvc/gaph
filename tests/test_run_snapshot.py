from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
ANALYSIS_DIR = REPO_ROOT / "analysis"
if str(ANALYSIS_DIR) not in sys.path:
    sys.path.insert(0, str(ANALYSIS_DIR))

from compare_runs_on_csv import build_outputs
from run_snapshot import export_runs, hydrate_exported_run_data, load_exported_run_data


def test_export_runs_and_load_exported_data(tmp_path: Path) -> None:
    source_root = tmp_path / "runs"
    run_dir = source_root / "group_a" / "run_alpha"
    _make_run(
        run_dir,
        params={"run_id": "alpha", "gene_ids": [1, 2], "read_generation": {"read_len": 80, "step": 40}},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIGCONF=Benign(1);AF=0.001;MC=SO:0001583|missense_variant;CSQ=missense_variant;HGVSC=c.100A>G",
                "1\t101\t.\tC\tT\t.\tPASS\tAF=0.2;CSQ=synonymous_variant",
            ],
        },
        failed_genes={"2": "No orthologs found"},
    )

    output_root = tmp_path / "exports"
    summaries = export_runs(source_root, output_root)

    assert len(summaries) == 1
    assert summaries[0].successful_genes == 1
    assert summaries[0].failed_genes == 1
    assert summaries[0].total_variants == 2

    export_dir = output_root / "group_a" / "run_alpha"
    assert (export_dir / "run_manifest.json").exists()
    assert (export_dir / "genes.csv.gz").exists()
    assert (export_dir / "variants.csv.gz").exists()
    assert (export_dir / "failure_events.csv.gz").exists()

    data = load_exported_run_data(export_dir, run_label="run_alpha")
    assert data["run_id"] == "run_alpha"
    assert len(data["gene_results"]) == 1
    assert data["gene_results"][0]["gene_id"] == "1"
    assert data["gene_results"][0]["total_variants"] == 2
    assert data["gene_results"][0]["clinvar_intersected"] == 1
    assert data["gene_results"][0]["gnomad_intersected"] == 2
    assert data["all_clnsig"] == ["Benign"]
    assert data["all_gnomad_af"] == [0.001, 0.2]

    manifest = json.loads((export_dir / "run_manifest.json").read_text())
    assert manifest["counts"]["successful_genes"] == 1
    assert manifest["counts"]["failed_genes"] == 1


def test_load_exported_run_data_summary_then_hydrate(tmp_path: Path) -> None:
    source_root = tmp_path / "runs"
    run_dir = source_root / "group_a" / "run_alpha"
    _make_run(
        run_dir,
        params={"run_id": "alpha", "gene_ids": [1]},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign;AF=0.001;CSQ=missense_variant",
            ],
        },
    )

    output_root = tmp_path / "exports"
    export_runs(source_root, output_root)
    export_dir = output_root / "group_a" / "run_alpha"

    summary = load_exported_run_data(export_dir, run_label="run_alpha", include_variants=False)
    assert summary["variant_rows"] == []
    assert summary["all_clnsig"] == []
    assert summary["gene_results"][0]["total_variants"] == 1
    assert summary["gene_results"][0]["clnsig_values"] == []

    hydrated = hydrate_exported_run_data(summary, export_dir)
    assert len(hydrated["variant_rows"]) == 1
    assert hydrated["all_clnsig"] == ["Benign"]
    assert hydrated["gene_results"][0]["clnsig_values"] == ["Benign"]


def test_export_treats_empty_final_vcf_as_incomplete(tmp_path: Path) -> None:
    source_root = tmp_path / "runs"
    run_dir = source_root / "run_alpha"
    _make_run(run_dir, params={"run_id": "alpha", "gene_ids": [1]})

    gene_dir = run_dir / "gene_1"
    gene_dir.mkdir(parents=True, exist_ok=True)
    (gene_dir / "gene_snps_annotated.vcf").write_text("")

    output_root = tmp_path / "exports"
    summaries = export_runs(source_root, output_root)

    assert summaries[0].successful_genes == 0
    assert summaries[0].incomplete_genes == 1

    data = load_exported_run_data(output_root / "run_alpha", run_label="run_alpha")
    assert data["gene_rows"][0]["status"] == "incomplete"
    assert data["gene_results"] == []


def test_compare_runs_on_exported_csv(tmp_path: Path) -> None:
    source_root = tmp_path / "runs"

    _make_run(
        source_root / "group_a" / "run_alpha",
        params={"run_id": "alpha", "gene_ids": [1], "read_generation": {"read_len": 80, "step": 40}},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIGCONF=Benign(1);AF=0.001;MC=SO:0001583|missense_variant;CSQ=missense_variant",
                "1\t101\t.\tC\tT\t.\tPASS\tAF=0.2;CSQ=synonymous_variant",
            ],
        },
    )
    _make_run(
        source_root / "group_a" / "run_beta",
        params={"run_id": "beta", "gene_ids": [1], "read_generation": {"read_len": 150, "step": 75}},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;AF=0.01;MC=SO:0001587|stop_gained;CSQ=stop_gained",
            ],
        },
    )

    output_root = tmp_path / "exports"
    export_runs(source_root, output_root)

    run_data = {
        "group_a/run_alpha": load_exported_run_data(output_root / "group_a" / "run_alpha", run_label="group_a/run_alpha"),
        "group_a/run_beta": load_exported_run_data(output_root / "group_a" / "run_beta", run_label="group_a/run_beta"),
    }
    results_dir = tmp_path / "results"
    outputs = build_outputs(run_data, results_dir)

    assert isinstance(outputs["df_stats"], pd.DataFrame)
    assert set(outputs["df_stats"]["Run"]) == {"group_a/run_alpha", "group_a/run_beta"}
    assert (results_dir / "overall_statistics.csv").exists()
    assert (results_dir / "params_comparison.csv").exists()
    assert (results_dir / "run_impact_class_shares.csv").exists()

    overall_stats = pd.read_csv(results_dir / "overall_statistics.csv")
    assert set(overall_stats["Run"]) == {"group_a/run_alpha", "group_a/run_beta"}


def test_compare_runs_on_precomputed_summaries(tmp_path: Path) -> None:
    source_root = tmp_path / "runs"

    _make_run(
        source_root / "group_a" / "run_alpha",
        params={"run_id": "alpha", "gene_ids": [1]},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIGCONF=Benign(1);AF=0.001;MC=SO:0001583|missense_variant;CSQ=missense_variant",
                "1\t101\t.\tC\tT\t.\tPASS\tAF=0.2;CSQ=synonymous_variant",
            ],
        },
    )
    _make_run(
        source_root / "group_a" / "run_beta",
        params={"run_id": "beta", "gene_ids": [1]},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;AF=0.01;MC=SO:0001587|stop_gained;CSQ=stop_gained",
            ],
        },
    )

    output_root = tmp_path / "exports"
    export_runs(source_root, output_root)

    run_data = {
        "group_a/run_alpha": load_exported_run_data(
            output_root / "group_a" / "run_alpha",
            run_label="group_a/run_alpha",
            include_variants=False,
        ),
        "group_a/run_beta": load_exported_run_data(
            output_root / "group_a" / "run_beta",
            run_label="group_a/run_beta",
            include_variants=False,
        ),
    }

    results_dir = tmp_path / "results_summary_only"
    outputs = build_outputs(run_data, results_dir)
    assert not outputs["df_clnsig"].empty
    assert not outputs["run_impact_class_shares"].empty
    assert (results_dir / "clnsig_distribution.csv").exists()


def _make_run(
    run_dir: Path,
    params: dict,
    genes: dict[str, list[str]] | None = None,
    failed_genes: dict[str, str] | None = None,
) -> None:
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "run_params.json").write_text(json.dumps(params))
    (run_dir / "ortholog_resolution.csv").write_text(
        "timestamp,gene_id,requested_scope,effective_scope,source_used,from_cache,species_count,sequence_count\n"
        "2026-01-01T00:00:00,1,all,all,ncbi_all,true,10,10\n"
    )

    genes = genes or {}
    for gene_id, variant_lines in genes.items():
        gene_dir = run_dir / f"gene_{gene_id}"
        gene_dir.mkdir(parents=True, exist_ok=True)
        header = [
            "##fileformat=VCFv4.3",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
        ]
        payload = "\n".join(header + variant_lines) + "\n"
        (gene_dir / "gene_snps_annotated.vcf").write_text(payload)

    failed_genes = failed_genes or {}
    if failed_genes:
        with (run_dir / "failed_genes.jsonl").open("w") as handle:
            for gene_id, message in failed_genes.items():
                gene_dir = run_dir / f"gene_{gene_id}"
                gene_dir.mkdir(parents=True, exist_ok=True)
                failure_payload = {
                    "gene_id": int(gene_id),
                    "gene_dir": str(gene_dir),
                    "status": "failed",
                    "started_at": "2026-01-01T00:00:00",
                    "finished_at": "2026-01-01T00:00:10",
                    "duration_seconds": 10.0,
                    "error_type": "RuntimeError",
                    "error_message": message,
                    "error_traceback": "Traceback",
                }
                handle.write(json.dumps(failure_payload) + "\n")
                (gene_dir / "failure.json").write_text(json.dumps(failure_payload))
