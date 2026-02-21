# Troubleshooting

Use this when genes fail or output looks suspicious.

## Failure Artifacts

Per failed gene:
- `runs/<run_id>/gene_<id>/failure.json` (`pipeline/pipeline.py:323`)

Run-level stream of failures:
- `runs/<run_id>/failed_genes.jsonl` (`pipeline/pipeline.py:324`)

Main run log:
- `runs/<run_id>/pipeline.log` (created at `pipeline/pipeline.py:227`)

Resume invocation log (if run was resumed):
- `runs/<run_id>/resume_invocations.jsonl`

## Common Failure Paths

### 1) No RefSeq proteins for BLAST fallback

Behavior:
- BLAST source requires `gene_protein_refseq` only (`pipeline/orthologs/blast_source.py:125`)
- Missing RefSeq IDs raises error (`pipeline/orthologs/blast_source.py:146`)

Effect:
- Gene fails and is recorded in `failure.json` / `failed_genes.jsonl`.

### 2) BAM filtering missing generated counts

Behavior:
- Generated counts are loaded from:
  - `genes_from_coordinates.fastq` first (`pipeline/bam_filtering.py:300`)
  - fallback `pseudo_reads.fastq` (`pipeline/bam_filtering.py:302`)
- If `min_mapped_pct_of_generated` is enabled and counts are missing, filtering raises runtime error (`pipeline/bam_filtering.py:506`).

### 3) gnomAD retrieval problems

Relevant code:
- NCBI->Ensembl candidates: `pipeline/gnomad.py:134`
- GraphQL fetch + retries: `pipeline/gnomad.py:314`

Check logs for selected Ensembl ID candidates and query errors.

### 4) Conda/mamba execution failures

All bio-tool commands go through `run_bio(...)` (`pipeline/config.py:121`).
On command failure, stderr/stdout are logged and exception is re-raised (`pipeline/config.py:127`).

### 5) Resume starts from unexpected gene

Behavior:
- Resume scanning is sequential over config `gene_ids` (`pipeline/pipeline.py:243`).
- Success is detected by `gene_snps_annotated.vcf`; failure by `failure.json` (`pipeline/pipeline.py:229`).
- Incomplete gene directories are deleted before rerun (`pipeline/pipeline.py:439`).

Checks:
1. Ensure `resume_run_dir` points to the intended run folder.
2. Verify config `gene_ids` order matches `run_params.json` from that run.
3. Inspect `resume_invocations.jsonl` and `pipeline.log` for selected `start_index` and `start_gene_id`.

## Debug Strategy

1. Read failure reason from `failure.json` first.
2. Jump to the referenced module/function using the pointers above.
3. Re-run one gene only (small config with one `gene_id`) to shorten feedback loop.
4. Verify final artifact contract remains intact: `gene_snps_annotated.vcf` exists on success.
