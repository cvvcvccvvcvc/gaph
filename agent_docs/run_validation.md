# Run and Validation

Use this for normal execution and quick verification after changes.

## Run

Standard run:

```bash
bash run.sh config.json
```

Resume an existing run:

```bash
# set "resume_run_dir" in config.json, then run the same command
bash run.sh config.json
```

Entrypoint details:
- env activation + python dispatch: `run.sh:31` to `run.sh:57`
- main function: `pipeline/pipeline.py:193`

## What to Verify (Minimum)

1. Run starts and writes run directory (`pipeline/pipeline.py:223`).
2. If `resume_run_dir` is set, verify logs show resume mode and start index (`pipeline/pipeline.py:358`).
3. Each successful gene ends with `gene_snps_annotated.vcf` (`pipeline/run_gene.py:561`).
4. If `keep_intermediate_files=false`, intermediates are removed (`pipeline/run_gene.py:565`).
5. No hidden gene failures: inspect `failed_genes.jsonl` if present (`pipeline/pipeline.py:279`).

## Useful Checks

Find latest run:

```bash
ls -lt runs | head
```

List failed genes quickly:

```bash
jq -r '.gene_id, .error_type, .error_message' runs/<run_id>/failed_genes.jsonl
```

Count successful final VCFs in run:

```bash
find runs/<run_id> -name gene_snps_annotated.vcf | wc -l
```

Inspect resume events:

```bash
cat runs/<run_id>/resume_invocations.jsonl
```

## Notes on Performance

- Ortholog prefetch is chunk-based in orchestrator: `pipeline/pipeline.py:257`, `pipeline/pipeline.py:281`
- Batch prefetch implementation is in `pipeline/orthologs/ncbi_datasets_source.py:207`
- Cache hit path for orthologs is in `pipeline/run_gene.py:394` and `pipeline/run_gene.py:407`
