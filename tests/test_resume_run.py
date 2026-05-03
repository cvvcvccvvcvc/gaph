from pathlib import Path
import csv
import gzip
import json
import sys
import types

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_DIR = REPO_ROOT / "pipeline"
if str(PIPELINE_DIR) not in sys.path:
    sys.path.insert(0, str(PIPELINE_DIR))

# Avoid importing heavy production dependencies for these unit tests.
bio_stub = types.ModuleType("Bio")
bio_stub.Entrez = types.SimpleNamespace(email=None, api_key=None)
sys.modules.setdefault("Bio", bio_stub)

run_gene_stub = types.ModuleType("run_gene")
run_gene_stub.run_gene = lambda *args, **kwargs: None
sys.modules.setdefault("run_gene", run_gene_stub)

orthologs_stub = types.ModuleType("orthologs")
orthologs_stub.get_source = lambda name: None
sys.modules.setdefault("orthologs", orthologs_stub)

class _LoggerStub:
    def remove(self, *args, **kwargs):
        return None

    def add(self, *args, **kwargs):
        return None

    def debug(self, *args, **kwargs):
        return None

    def info(self, *args, **kwargs):
        return None

    def warning(self, *args, **kwargs):
        return None

    def success(self, *args, **kwargs):
        return None

    def error(self, *args, **kwargs):
        return None

    def exception(self, *args, **kwargs):
        return None


loguru_stub = types.ModuleType("loguru")
loguru_stub.logger = _LoggerStub()
sys.modules.setdefault("loguru", loguru_stub)

import pipeline as pipeline_main  # noqa: E402
import run_compaction  # noqa: E402
from run_compaction import compact_run_in_place, load_exported_run_data  # noqa: E402


def _touch(path: Path, content: str = "x") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


def test_validate_resume_run_dir_cfg():
    assert pipeline_main._validate_resume_run_dir_cfg({}) is None
    assert pipeline_main._validate_resume_run_dir_cfg({"resume_run_dir": "runs/run_1"}) == "runs/run_1"

    with pytest.raises(TypeError):
        pipeline_main._validate_resume_run_dir_cfg({"resume_run_dir": ""})

    with pytest.raises(TypeError):
        pipeline_main._validate_resume_run_dir_cfg({"resume_run_dir": 123})


def test_validate_output_compaction_cfg():
    assert pipeline_main.validate_output_compaction_cfg({}) == {"enabled": False}
    assert pipeline_main.validate_output_compaction_cfg({"output_compaction": None}) == {"enabled": False}
    assert pipeline_main.validate_output_compaction_cfg({"output_compaction": {"enabled": True}}) == {"enabled": True}

    with pytest.raises(TypeError):
        pipeline_main.validate_output_compaction_cfg({"output_compaction": []})

    with pytest.raises(TypeError):
        pipeline_main.validate_output_compaction_cfg({"output_compaction": {"enabled": "yes"}})


def test_gene_status_success_failed_incomplete_missing(tmp_path):
    run_dir = tmp_path / "run_1"

    _touch(run_dir / "gene_1" / "gene_snps_annotated.vcf", "##fileformat=VCFv4.2\n")
    _touch(run_dir / "gene_2" / "failure.json", '{"status":"failed"}')
    _touch(run_dir / "gene_3" / "partial.tmp", "work in progress")

    assert pipeline_main._gene_status(run_dir, 1) == "success"
    assert pipeline_main._gene_status(run_dir, 2) == "failed"
    assert pipeline_main._gene_status(run_dir, 3) == "incomplete"
    assert pipeline_main._gene_status(run_dir, 4) == "missing"


def test_gene_status_uses_compacted_snapshot_when_gene_dir_was_deleted(tmp_path):
    run_dir = tmp_path / "run_1"
    run_dir.mkdir()
    (run_dir / "run_manifest.json").write_text(
        json.dumps({"counts": {"gene_rows": 2, "successful_genes": 1, "failed_genes": 1}})
    )
    with gzip.open(run_dir / "genes.csv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["gene_id", "status"])
        writer.writeheader()
        writer.writerow({"gene_id": "1", "status": "success"})
        writer.writerow({"gene_id": "2", "status": "failed"})

    assert pipeline_main._gene_status(run_dir, 1) == "success"
    assert pipeline_main._gene_status(run_dir, 2) == "failed"
    assert pipeline_main._gene_status(run_dir, 3) == "missing"


def test_compact_run_in_place_writes_snapshot_and_deletes_gene_dirs(tmp_path):
    run_dir = tmp_path / "run_alpha"
    _make_compaction_run(
        run_dir,
        params={"run_id": "alpha", "gene_ids": [1, 2]},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tCLNSIG=Benign;AF=0.001;CSQ=missense_variant",
            ],
        },
        failed_genes={"2": "No orthologs found"},
    )
    stale_temp_dir = run_dir / ".compaction_stale"
    stale_temp_dir.mkdir()
    (stale_temp_dir / "partial.tmp").write_text("left by interrupted compaction")

    result = compact_run_in_place(run_dir)

    assert result.successful_genes == 1
    assert result.failed_genes == 1
    assert result.deleted_gene_dirs == 2
    assert not (run_dir / "gene_1").exists()
    assert not (run_dir / "gene_2").exists()
    assert not stale_temp_dir.exists()
    assert (run_dir / "genes.csv.gz").exists()
    assert (run_dir / "variants.csv.gz").exists()
    assert (run_dir / "failure_events.csv.gz").exists()
    assert (run_dir / "analysis_summary.json.gz").exists()

    manifest = json.loads((run_dir / "run_manifest.json").read_text())
    assert manifest["compaction"]["raw_gene_dirs_deleted"] is True
    assert manifest["counts"]["successful_genes"] == 1
    assert manifest["counts"]["failed_genes"] == 1

    with gzip.open(run_dir / "analysis_summary.json.gz", "rt") as handle:
        summary = json.load(handle)
    assert summary["clnsig_counts"] == [{"CLNSIG": "Benign", "count": 1}]
    assert summary["af_stats"]["count"] == 1
    assert summary["run_impact_rows"] == [
        {
            "source": "gnomAD_CSQ",
            "impact_class": "MODERATE",
            "count": 1,
            "source_total": 1,
            "share_%": 100.0,
        }
    ]

    data = load_exported_run_data(run_dir, run_label="run_alpha")
    assert len(data["gene_rows"]) == 2
    assert len(data["variant_rows"]) == 1


def test_compact_run_in_place_does_not_load_variants_during_validation(tmp_path, monkeypatch):
    run_dir = tmp_path / "run_alpha"
    _make_compaction_run(
        run_dir,
        params={"run_id": "alpha", "gene_ids": [1]},
        genes={
            "1": [
                "1\t100\t.\tA\tG\t.\tPASS\tAF=0.001;CSQ=missense_variant",
                "1\t101\t.\tC\tT\t.\tPASS\tAF=0.2;CSQ=synonymous_variant",
            ],
        },
    )

    original_read_gzip_csv = run_compaction._read_gzip_csv

    def guarded_read_gzip_csv(path: Path):
        if Path(path).name == "variants.csv.gz":
            raise AssertionError("variants.csv.gz should be counted as a stream during validation")
        return original_read_gzip_csv(path)

    monkeypatch.setattr(run_compaction, "_read_gzip_csv", guarded_read_gzip_csv)

    result = run_compaction.compact_run_in_place(run_dir)

    assert result.total_variants == 2
    assert not (run_dir / "gene_1").exists()


def test_compact_run_in_place_refuses_to_delete_incomplete_runs(tmp_path):
    run_dir = tmp_path / "run_alpha"
    _make_compaction_run(
        run_dir,
        params={"run_id": "alpha", "gene_ids": [1, 2]},
        genes={"1": ["1\t100\t.\tA\tG\t.\tPASS\tAF=0.001"]},
    )
    incomplete_dir = run_dir / "gene_2"
    incomplete_dir.mkdir(parents=True)
    (incomplete_dir / "partial.tmp").write_text("in progress")

    with pytest.raises(ValueError, match="Refusing to delete raw outputs"):
        compact_run_in_place(run_dir)

    assert (run_dir / "gene_1" / "gene_snps_annotated.vcf").exists()
    assert (run_dir / "gene_2" / "partial.tmp").exists()
    assert not (run_dir / "run_manifest.json").exists()


def test_build_resume_plan_starts_from_first_incomplete(tmp_path):
    run_dir = tmp_path / "run_1"

    _touch(run_dir / "gene_10" / "gene_snps_annotated.vcf", "done")
    _touch(run_dir / "gene_20" / "pseudo_reads.fastq", "partial")

    plan = pipeline_main._build_resume_plan(run_dir, [10, 20, 30, 40])

    assert plan["start_index"] == 1
    assert plan["start_status"] == "incomplete"
    assert plan["counts"] == {
        "success": 1,
        "failed": 0,
        "incomplete": 1,
        "missing": 2,
    }


def test_build_resume_plan_starts_from_first_missing_when_prior_are_terminal(tmp_path):
    run_dir = tmp_path / "run_1"

    _touch(run_dir / "gene_10" / "gene_snps_annotated.vcf", "done")
    _touch(run_dir / "gene_20" / "failure.json", '{"status":"failed"}')

    plan = pipeline_main._build_resume_plan(run_dir, [10, 20, 30, 40])

    assert plan["start_index"] == 2
    assert plan["start_status"] == "missing"
    assert plan["counts"] == {
        "success": 1,
        "failed": 1,
        "incomplete": 0,
        "missing": 2,
    }


def test_build_resume_plan_returns_none_if_all_genes_terminal(tmp_path):
    run_dir = tmp_path / "run_1"

    _touch(run_dir / "gene_10" / "gene_snps_annotated.vcf", "done")
    _touch(run_dir / "gene_20" / "failure.json", '{"status":"failed"}')

    plan = pipeline_main._build_resume_plan(run_dir, [10, 20])

    assert plan["start_index"] is None
    assert plan["start_status"] is None
    assert plan["counts"] == {
        "success": 1,
        "failed": 1,
        "incomplete": 0,
        "missing": 0,
    }


def test_read_run_params_gene_ids(tmp_path):
    run_dir = tmp_path / "run_1"
    _touch(run_dir / "run_params.json", '{"gene_ids": [1, "2", 3]}')
    assert pipeline_main._read_run_params_gene_ids(run_dir) == [1, 2, 3]

    _touch(run_dir / "run_params.json", '{"gene_ids": "not-a-list"}')
    with pytest.raises(TypeError):
        pipeline_main._read_run_params_gene_ids(run_dir)


def _make_compaction_run(
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
