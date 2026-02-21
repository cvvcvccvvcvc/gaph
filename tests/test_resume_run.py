from pathlib import Path
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


def test_gene_status_success_failed_incomplete_missing(tmp_path):
    run_dir = tmp_path / "run_1"

    _touch(run_dir / "gene_1" / "gene_snps_annotated.vcf", "##fileformat=VCFv4.2\n")
    _touch(run_dir / "gene_2" / "failure.json", '{"status":"failed"}')
    _touch(run_dir / "gene_3" / "partial.tmp", "work in progress")

    assert pipeline_main._gene_status(run_dir, 1) == "success"
    assert pipeline_main._gene_status(run_dir, 2) == "failed"
    assert pipeline_main._gene_status(run_dir, 3) == "incomplete"
    assert pipeline_main._gene_status(run_dir, 4) == "missing"


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
