# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a bioinformatics pipeline for analyzing gene variants using evolutionary conservation. The project uses BRCA1 (Gene ID 672) as the primary example but the pipeline is designed to be generalizable to other genes.

The workflow identifies potentially pathogenic variants by:
1. Finding homologous gene sequences across species (via BLAST or NCBI Datasets CLI)
2. Aligning homolog sequences to the human reference
3. Calling variants at conserved/divergent positions
4. Annotating variants with ClinVar and gnomAD databases to validate the approach

## Environment Setup

```bash
# Single conda env with Python + bio tools (see environment.yml)
micromamba create -f environment.yml   # or: conda env create -f environment.yml

# Run the pipeline
./run.sh config.json

# Or manually:
micromamba activate bio
python pipeline/pipeline.py config.json
```

All dependencies (Python packages + bioinformatics tools + NCBI datasets CLI) are in the `bio` conda environment defined by `environment.yml`.

### Configuration

All settings are in `config.json`:

```json
{
  "gene_ids": [672],
  "hitlist_size": 5000,
  "env_file": "~/.secrets.env",
  "data_dir": "data",
  "runs_dir": "runs",
  "gnomad_dir": "gnomad_variants",
  "conda_env": "bio",
  "keep_intermediate_files": false,
  "bam_filtering": {
    "enabled": true,
    "min_mapped_pct_of_generated": 10,
    "max_pct_filtered": 50,
    "min_kept_pct_of_reference": 10,
    "read_len": 75,
    "step": 35
  }
}
```

- `env_file` — path to `.env` file containing `ENTREZ_EMAIL` and `ENTREZ_API_KEY` (tilde-expanded). Secrets stay out of config.json.
- `data_dir`, `runs_dir`, `gnomad_dir` — relative to repo root, or absolute paths
- `conda_env` — name of the conda environment (default: `bio`)
- `keep_intermediate_files` — keep all per-gene intermediates when `true`; default `false` keeps only `gene_snps_annotated.vcf` in each gene output folder
- `bam_filtering` — LIS-based read-order filtering before variant calling
  - `enabled` — enable/disable filtering
  - `min_mapped_pct_of_generated`, `max_pct_filtered`, `min_kept_pct_of_reference` — required when enabled, each in `[0,100]`
  - `read_len`, `step` — pseudo-read generation parameters for denominator calculations (default `75` and `35`)

If `bam_filtering.enabled=true` and `min_mapped_pct_of_generated` is set, missing generated pseudo-read counts for any homologue cause a fail-fast error for that gene.

One-time cluster setup: `./setup_cluster.sh` (creates dirs, downloads ClinVar)

## Key Commands

### BWA Alignment
```bash
micromamba run -n bio bwa index reference.fasta
micromamba run -n bio bwa mem reference.fasta reads.fastq > aln.sam
micromamba run -n bio samtools view -bS aln.sam > aln.bam
micromamba run -n bio samtools sort aln.bam -o aln.sorted.bam
micromamba run -n bio samtools index aln.sorted.bam
```

### Variant Calling (VarScan)
```bash
micromamba run -n bio samtools faidx reference.fasta
micromamba run -n bio samtools mpileup -f reference.fasta aln.sorted.bam > pileup.mpileup
micromamba run -n bio varscan mpileup2snp pileup.mpileup --min-coverage 8 --min-reads2 2 --output-vcf 1 > snps.vcf
```

### VCF Processing & Annotation
```bash
# Compress and index normalized VCF
micromamba run -n bio bgzip -c normalized.vcf > normalized.vcf.gz
micromamba run -n bio tabix -p vcf normalized.vcf.gz

# Annotate with ClinVar and gnomAD in a single pass
micromamba run -n bio bcftools annotate -a clinvar.vcf.gz -c INFO/CLNSIG,INFO/CLNDN,INFO/CLNREVSTAT sample.vcf.gz \
  | micromamba run -n bio bcftools annotate -a gnomad.vcf.gz -c INFO/AF,INFO/MAF,INFO/CSQ -o annotated.vcf -O v
```

## Project Architecture

### Directory Structure
- `data/` - Shared reference data (ClinVar VCF)
- `gnomad_variants/` - Cached gnomAD VCF files by Ensembl ID
- `runs/` - Per-run outputs (was `pipeline/runs/`)
  - `runs/run_YYYYMMDD_HHMMSS/` - Timestamped run folder
    - `run_params.json` - Run configuration parameters
    - `pipeline.log` - Run log file
    - `gene_<ID>/` - Per-gene outputs
      - `gene_seq.fasta` - Reference gene sequence
      - `genes_from_coordinates.fastq` - Homolog sequences
      - `pseudo_reads.fastq` - Generated pseudoreads
      - `aln.sorted.bam` - Aligned reads
      - `aln.filtered.lis.bam` - Filtered reads (when BAM filtering enabled and intermediates kept)
      - `bam_filtering_stats.json` - Per-homologue BAM filtering stats (when intermediates kept)
      - `bam_filtering_overall.json` - BAM filtering summary (when intermediates kept)
      - `gene_snps.vcf` - Called variants (local coordinates)
      - `gene_snps_normalized.vcf.gz` - Normalized variants (genomic coordinates)
      - `gene_snps_annotated.vcf` - Variants annotated with ClinVar and gnomAD
      - `ncbi_dataset/` - NCBI Datasets CLI output (if using ncbi_datasets source)
- `pipeline/` - Main pipeline code and modules
  - `pipeline/config.py` - Central configuration (paths, conda, credentials)
  - `pipeline/pipeline.py` - Orchestrator: reads config, loops genes
  - `pipeline/run_gene.py` - Process one gene through all pipeline steps
  - `pipeline/gnomad.py` - gnomAD integration module
  - `pipeline/orthologs/` - Ortholog sources (BLAST, NCBI Datasets CLI)
- `_cluster/` - Local cluster simulation for testing
- `run.sh` - Bash entry point (read config.json, activate env, run pipeline)
- `setup_cluster.sh` - One-time cluster setup (create dirs, download ClinVar)
- `environment.yml` - Conda environment definition
- `config.json.example` - Pipeline config template
- `analyze_results.ipynb` - Cross-gene analysis and visualization

### Main Files
- `pipeline/pipeline.py` - Main orchestrator (entry point)
- `pipeline/run_gene.py` - All per-gene pipeline steps
- `pipeline/config.py` - Central configuration
- `pipeline/gnomad.py` - gnomAD integration module
- `pipeline/orthologs/` - Pluggable ortholog source system
- `analyze_results.ipynb` - Cross-run analysis for comparing variants across genes

### Data Flow
1. **Gene ID → Genomic Sequence**: Fetch human gene sequence and coordinates from NCBI
2. **Gene ID → Orthologs**: Find homologs via BLAST or NCBI Datasets CLI
3. **Orthologs → Genomic Sequences → FASTQ**: Fetch homolog gene sequences from NCBI
4. **FASTQ → Pseudoreads**: Slice sequences into overlapping reads (75bp, step 35)
5. **Pseudoreads → Alignment → BAM**: Align homolog reads to human reference with BWA
6. **BAM → Filtered BAM**: Filter by dominant strand + LIS backbone + homologue-level thresholds
7. **Filtered BAM → Variants → VCF**: Call SNPs/indels with VarScan
8. **VCF → Normalized VCF**: Convert local coordinates to genomic coordinates
9. **Normalized VCF + ClinVar + gnomAD → Annotated VCF**: Single-pass annotation with both databases

### Key Functions (in pipeline/run_gene.py)

All functions take an explicit `work_dir: Path` parameter — no `os.chdir()`.

**Core Pipeline Functions:**
- `download_gene_seq(gene_id, work_dir)` - Fetch genomic sequence, returns coordinates dict
- `generate_pseudoreads(input_fastq, work_dir)` - Create sliding window reads from sequences (75bp, step 35)
- `align_pseudoreads(work_dir)` - Run BWA alignment pipeline (index, align, sort, index BAM)
- `filter_bam_for_gene(work_dir, filtering_cfg, verbose=False)` - Production LIS BAM filtering
- `call_variants(work_dir, bam_path)` - Run VarScan variant calling (SNPs and indels)
- `normalize_vcf(gene_coords, work_dir)` - Convert local VCF coordinates to genomic coordinates
- `annotate_variants(work_dir, gnomad_vcf_gz)` - Annotate variants with ClinVar and gnomAD in a single pass
- `run_gene(gene_id, work_dir, cfg)` - Execute full pipeline for one gene

**Ortholog System (pipeline/orthologs/):**
- `get_source(name)` - Get ortholog source by name ("blast" or "ncbi_datasets")
- `list_sources()` - List available ortholog sources
- `source.fetch(gene_id, output_dir)` - Fetch orthologs, returns OrthologResult with fastq_path

**gnomAD Module (pipeline/gnomad.py):**
- `ncbi_to_ensembl(ncbi_gene_id)` - Convert NCBI Gene ID to Ensembl ID via MyGene.info
- `download_gnomad_for_gene(ncbi_gene_id)` - Download gnomAD variants (cached), returns VCF path
- `prepare_gnomad_vcf(vcf_path, output_dir)` - Compress and index gnomAD VCF for bcftools

## Analysis & Visualization

The [analyze_results.ipynb](analyze_results.ipynb) notebook provides cross-gene analysis:
- Collects variant counts from all processed genes in a run
- Calculates ClinVar and gnomAD intersection rates
- Parses and visualizes clinical significance (CLNSIG) distribution
- Analyzes gnomAD allele frequency distributions
- Creates summary tables and plots for publication
- Exports results to CSV files in `analysis_results/`

## Notes

- NCBI API requires email and API key (set in Entrez.email and Entrez.api_key)
- ClinVar file is shared across all genes: `data/clinvar.vcf.gz`
- gnomAD files are cached per-gene in `gnomad_variants/` by Ensembl ID
- gnomAD download uses MyGene.info API for NCBI→Ensembl ID conversion
- gnomAD GraphQL API has rate limits (6 second delay between queries)
- VCF normalization converts region-relative coordinates to chromosome coordinates
- Some genes lack genomic coordinates in NCBI - these are skipped during processing
- Logging via loguru outputs to both stdout and per-run log file inside `runs/run_*/pipeline.log`
- Each pipeline run writes `run_params.json` in the run folder to capture parameters
- Ortholog sources are pluggable - use `get_source("blast")` or `get_source("ncbi_datasets")`
