# Project Map

This file is a quick code-navigation map for AI agents.

## Production Boundary

Use `pipeline/` for production logic.
Use `tests/` for experiments, plots, and sandbox scripts.
Do not move test-only helpers into `pipeline/` unless explicitly requested.

## Entrypoints

- Shell entrypoint: `run.sh:1`
- Python run orchestrator: `pipeline/pipeline.py:193`
- Per-gene execution: `pipeline/run_gene.py:494`

## Core Modules

- Runtime config/env + bio command runner: `pipeline/config.py:81`, `pipeline/config.py:121`
- Resume planning and gene terminal-state detection: `pipeline/pipeline.py:229`, `pipeline/pipeline.py:243`
- Ortholog fetch (cache-first + source fallback): `pipeline/run_gene.py:390`
- NCBI Datasets ortholog source + batch prefetch: `pipeline/orthologs/ncbi_datasets_source.py:169`, `pipeline/orthologs/ncbi_datasets_source.py:207`
- BLAST ortholog source (RefSeq proteins only): `pipeline/orthologs/blast_source.py:37`, `pipeline/orthologs/blast_source.py:123`
- BAM filtering (LIS production path): `pipeline/bam_filtering.py:353`
- gnomAD fetch + NCBI->Ensembl candidate logic: `pipeline/gnomad.py:134`, `pipeline/gnomad.py:314`

## Output Contract

- Final expected per-gene artifact: `gene_snps_annotated.vcf` (`pipeline/run_gene.py:561`)
- If `keep_intermediate_files=false`, gene folder is cleaned to final output only (`pipeline/run_gene.py:565`, `pipeline/run_gene.py:473`)
- Failed genes are recorded as:
  - per-gene: `failure.json` (`pipeline/pipeline.py:323`)
  - run-level append log: `failed_genes.jsonl` (`pipeline/pipeline.py:279`, `pipeline/pipeline.py:324`)
