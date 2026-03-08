# Bioinf Variant Pipeline

Pipeline for gene-level variant discovery from ortholog-derived pseudo-reads, with ClinVar and gnomAD annotation.

## What it does

1. Downloads reference gene sequence (NCBI).
2. Fetches ortholog sequences (NCBI Datasets, fallback to BLAST).
3. Generates pseudo-reads from ortholog sequences.
4. Aligns pseudo-reads to reference gene (BWA + samtools).
5. Filters BAM by homologue read-order consistency (LIS-based filtering).
6. Calls variants (VarScan).
7. Normalizes VCF coordinates to genomic coordinates.
8. Annotates variants with ClinVar and optionally gnomAD.

## Requirements

- micromamba/mamba/conda
- environment from `environment.yml`
- `config.json` in repo root

## Setup

```bash
micromamba create -f environment.yml
```

## Run

```bash
bash run.sh config.json
```

## Configuration (`config.json`)

Example:

```json
{
  "gene_ids": [672],
  "hitlist_size": 5000,
  "blast_expect": 10.0,
  "read_generation": {
    "pseudo_read_phred": 30,
    "read_len": 75,
    "step": 35
  },
  "variant_calling": {
    "min_var_freq": 0.2,
    "min_coverage": 8,
    "min_reads2": 2
  },
  "ortholog_selection": {
    "scope": "all"
  },
  "env_file": ".env",
  "data_dir": "data",
  "runs_dir": "runs",
  "resume_run_dir": "runs/run_20260221_125031",
  "gnomad_dir": "gnomad_variants",
  "conda_env": "bio",
  "keep_intermediate_files": false,
  "bam_filtering": {
    "enabled": true,
    "wrong_strand": true,
    "lis": true,
    "overlap": true,
    "min_mapped_pct_of_generated": 10,
    "max_pct_filtered": 50,
    "min_kept_pct_of_reference": 10
  }
}
```

### Resume an existing run

- `resume_run_dir` is optional.
- If omitted, pipeline creates a new `runs/run_YYYYMMDD_HHMMSS` directory.
- If provided, pipeline inspects configured genes in order and resumes from:
  - the first gene with a non-terminal state (`gene_<id>/` exists but no `gene_snps_annotated.vcf` and no `failure.json`), or
  - the first missing gene directory if all previous genes are terminal.
- Terminal state is one of:
  - `gene_snps_annotated.vcf` exists and is non-empty (success)
  - `failure.json` exists and is non-empty (failed)
- If all genes are terminal in `resume_run_dir`, a new run directory is created automatically.
- In resume mode, existing terminal genes are skipped; incomplete gene directories are deleted and recomputed from scratch.

### `bam_filtering` rules

- `enabled`: boolean.
- `wrong_strand`: boolean (enable/disable dominant-strand filtering stage).
- `lis`: boolean (enable/disable LIS/LDS ordering stage).
- `overlap`: boolean (enable/disable overlap dedup stage).
- If `enabled=true`, these are required and must be in `[0, 100]`:
  - `min_mapped_pct_of_generated`
  - `max_pct_filtered`
  - `min_kept_pct_of_reference`
- If `min_mapped_pct_of_generated` is enabled and generated pseudo-read counts cannot be computed for homologues, gene processing fails fast.

### BLAST parameters

- `hitlist_size`: max number of BLAST hits to request.
- `blast_expect`: BLAST E-value cutoff (passed to `qblast(..., expect=...)`).

### Ortholog selection

- `ortholog_selection.scope`: requested NCBI ortholog scope (`all`, `mammals`, `eukaryota`, etc.).
- Resolution order is fixed:
  1. NCBI Datasets with requested scope
  2. NCBI Datasets with `all`
  3. BLAST fallback
- NCBI ortholog cache is scope-aware when `scope != all`.

### Read generation parameters

- `read_generation.pseudo_read_phred`: fixed PHRED assigned to generated pseudo-reads.
- `read_generation.read_len`: pseudo-read length.
- `read_generation.step`: sliding-window step for pseudo-read generation.

### Variant calling parameters

- `variant_calling.min_var_freq`: passed to VarScan as `--min-var-freq`.
- `variant_calling.min_coverage`: passed to VarScan as `--min-coverage`.
- `variant_calling.min_reads2`: passed to VarScan as `--min-reads2`.

## Outputs

Run folder:

- `runs/run_YYYYMMDD_HHMMSS/pipeline.log`
- `runs/run_YYYYMMDD_HHMMSS/run_params.json`
- `runs/run_YYYYMMDD_HHMMSS/ortholog_resolution.csv`
- `runs/run_YYYYMMDD_HHMMSS/gene_<ID>/gene_snps_annotated.vcf` (final output)

If `keep_intermediate_files=true`, gene directory also keeps intermediates such as:

- `aln.sorted.bam`
- `aln.filtered.lis.bam`
- `bam_filtering_stats.json`
- `bam_filtering_overall.json`
- other intermediate FASTQ/VCF/BAM files

## Notes

- Entrez credentials are loaded from `env_file` (`ENTREZ_EMAIL`, `ENTREZ_API_KEY`).
- gnomAD files are cached in `gnomad_variants/`.
- Test/sandbox plotting tools stay under `tests/` and are not part of production pipeline execution.
