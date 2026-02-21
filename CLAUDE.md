# CLAUDE.md

This file is the minimal always-loaded onboarding for AI agents in this repo.

## WHY

Goal of the project: run a per-gene variant pipeline from ortholog-derived pseudo-reads and produce `gene_snps_annotated.vcf` for each gene.

Primary production path:
1. ortholog retrieval
2. pseudoread generation
3. BWA/SAMtools alignment
4. optional BAM filtering
5. VarScan calling
6. ClinVar + gnomAD annotation

## WHAT

Core code:
- `pipeline/pipeline.py` - run orchestration, config validation, run-level logging.
- `pipeline/run_gene.py` - per-gene pipeline steps and cleanup contract.
- `pipeline/config.py` - path/env initialization and `run_bio(...)`.
- `pipeline/bam_filtering.py` - production LIS-based BAM filtering.
- `pipeline/orthologs/` - ortholog sources (`ncbi_datasets`, `blast`).
- `pipeline/gnomad.py` - NCBI->Ensembl mapping and gnomAD fetching.

Operational entrypoint:
- `run.sh` -> activates conda/mamba env and runs `pipeline/pipeline.py`.

## HOW

Default execution:
- `bash run.sh config.json`

Agent workflow rules:
1. Keep this file short and universal. Put task-specific detail in `agent_docs/`.
2. Before coding, read only relevant docs from `agent_docs/` (do not load all by default).
3. Prefer pointers to authoritative code (`file:line`) over duplicated snippets.
4. Do not add style-guide/lint rules here; rely on project code patterns and normal tooling.
5. Never commit secrets. Credentials come from `env_file` configured in `config.json`.
6. Preserve production vs sandbox split: production logic lives in `pipeline/`, experiments/plots in `tests/`.
7. For large new features, create and use a dedicated git branch (do not implement them directly on `main`).
8. Keep commits granular: commit after each completed fix/step with a focused diff.

## Progressive Disclosure

Read these only when relevant:

- `agent_docs/project_map.md` - architecture map and file ownership.
- `agent_docs/config_contract.md` - config schema, validation rules, runtime contracts.
- `agent_docs/run_validation.md` - run commands and practical verification flow.
- `agent_docs/troubleshooting.md` - failure artifacts and common incident paths.
