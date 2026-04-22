# Config Contract

Use this when changing `config.json` schema or behavior.

## Where Config Is Validated

Validation happens in `pipeline/pipeline.py` before processing genes:
- numeric core fields: hitlist/blast parameters in `pipeline/pipeline.py`
- `read_generation` block: `pipeline/pipeline.py`
- `variant_calling` block: `pipeline/pipeline.py`
- `ortholog_selection` block: `pipeline/pipeline.py`
- resume run path: `pipeline/pipeline.py:197`, `pipeline/pipeline.py:296`
- BAM filtering block: `pipeline/pipeline.py:86`
- cache block: `pipeline/pipeline.py:135`
- output compaction block: `pipeline/pipeline.py`

If you add a config field, update validation here first.

## Runtime Initialization

Resolved paths, env-file loading, and conda resolver are in `pipeline/config.py:81`.
Critical details:
- `env_file` is parsed in `pipeline/config.py:93`
- project paths are resolved relative to repo root via `pipeline/config.py:33`
- bio commands run through `run_bio(...)` in `pipeline/config.py:121`
- `resume_run_dir` (if set) is resolved relative to repo root via `pipeline/pipeline.py:206`
- `run_name` (if set) is validated in `pipeline/pipeline.py` and used as the readable prefix for new run directories

## Resume Rules

- Optional config field: `resume_run_dir` (string path to existing run folder).
- Optional config field: `run_name` (string prefix for new run folder names).
- Gene terminal states in resume mode (`pipeline/pipeline.py:229`):
  - success: `gene_snps_annotated.vcf` exists and non-empty
  - failed: `failure.json` exists and non-empty
  - compacted success/failed: status row in `genes.csv.gz` when raw gene files were deleted
- Resume start point is first non-terminal gene in config order (`pipeline/pipeline.py:243`).
- If all genes are terminal, pipeline starts a new run directory (`pipeline/pipeline.py:316`).
- A compacted run with incomplete/missing genes is not resumable.
- New run directories keep the `run_` prefix and use `run_name` or the config stem as the base name.
- If the base name already exists, a numeric suffix is appended: `__2`, `__3`, ...

## Output Compaction Rules

- Optional config block: `output_compaction`.
- `output_compaction.enabled` defaults to `false`.
- When enabled, compaction runs after all genes finish and writes compact run-level files:
  - `run_manifest.json`
  - `genes.csv.gz`
  - `variants.csv.gz`
  - `failure_events.csv.gz`
  - `analysis_summary.json.gz`
- Raw `gene_<ID>/` directories are deleted only after compact files are validated.

## BAM Filtering Rules

- `bam_filtering` block is required by current validator (`pipeline/pipeline.py:88`)
- per-stage toggles: `wrong_strand`, `lis`, `overlap` (booleans)
- required thresholds when enabled: `pipeline/pipeline.py:108`
- pseudo-read geometry (`read_len`, `step`) now lives in `read_generation` and is passed into filtering runtime from `pipeline/run_gene.py`
- fail-fast on missing generated counts is enforced in filtering module: `pipeline/bam_filtering.py:506`

## Cache Rules

- cache flags are normalized in `pipeline/pipeline.py:159`
- cache paths are expanded in `pipeline/config.py:108` to `pipeline/config.py:114`
- NCBI ortholog cache keys are scope-aware (`all` keeps legacy file names, other scopes use `__scope_<slug>`)
- ortholog cache writes currently use gzip `compresslevel=1`:
  - `pipeline/run_gene.py:86`
  - `pipeline/orthologs/ncbi_datasets_source.py:121`

## Ortholog Resolution Rules

- Requested scope comes from `ortholog_selection.scope`.
- Resolution chain in runtime is:
  1. NCBI requested scope (cache, then live fetch)
  2. NCBI `all` (cache, then live fetch)
  3. BLAST (cache, then live fetch)
- Run-level audit is written to `ortholog_resolution.csv` in the run directory.

## Change Checklist

1. Update validator(s) in `pipeline/pipeline.py`.
2. Update runtime usage in `pipeline/config.py` and/or `pipeline/run_gene.py`.
3. Keep `config.json` example(s) aligned.
4. Run a small real pipeline run and verify `run_params.json` captures new fields (`pipeline/pipeline.py:230`).
