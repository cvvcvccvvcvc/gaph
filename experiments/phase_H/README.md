# Phase I: VarScan threshold sweep

Goal: tune VarScan after fixing the upstream pipeline at the current balanced
candidate:

- read geometry: `300/30`
- ortholog scope: `mammals`
- BAM filtering: `wrong_strand + LIS`, `overlap=false`
- BAM thresholds: `min_mapped_pct_of_generated=0`, `max_pct_filtered=100`,
  `min_kept_pct_of_reference=0`

`run_36_min_mapped_pct_0` is the external control for this phase:
`min_coverage=1`, `min_reads2=1`, `min_var_freq=0.01`.

## Why this phase

Previous experiments already show the main upstream shape:

- `300/30` captures almost all useful read-geometry signal at much lower cost
  than `500/30`.
- `mammals` is much faster than `all` and has the best benign/pathogenic
  trade-off among tested scopes.
- LIS is the useful BAM filter; `wrong_strand` is kept because it is cheap and
  does not materially change the callset after LIS.
- `min_mapped_pct_of_generated` above `10%` starts removing too much ClinVar
  signal, so VarScan should be tuned before adding stricter BAM thresholds.

The sweep is sparse by design. It separates the early effects of coverage,
alternate-read support, and variant allele fraction before moving to stricter
combined thresholds.

## Configs

| config | min_coverage | min_reads2 | min_var_freq | purpose |
| --- | ---: | ---: | ---: | --- |
| `37_varscan_cov2_reads1_vaf001.json` | 2 | 1 | 0.01 | first coverage gate above recall control |
| `38_varscan_cov4_reads1_vaf001.json` | 4 | 1 | 0.01 | coverage-only pressure without read-support pressure |
| `39_varscan_cov4_reads2_vaf001.json` | 4 | 2 | 0.01 | isolate `min_reads2` at same coverage and VAF |
| `40_varscan_cov4_reads2_vaf002.json` | 4 | 2 | 0.02 | very mild VAF gate |
| `41_varscan_cov4_reads2_vaf005.json` | 4 | 2 | 0.05 | balanced early production candidate |
| `42_varscan_cov6_reads2_vaf005.json` | 6 | 2 | 0.05 | coverage-tightened version of 41 |
| `43_varscan_cov8_reads2_vaf010.json` | 8 | 2 | 0.10 | conservative midpoint |
| `44_varscan_cov8_reads2_vaf020.json` | 8 | 2 | 0.20 | default-like strict VarScan anchor |

Recommended first batch if compute is limited: run `37` through `41`, compare
against `run_36`, then run `42` through `44` only if the curve still has no
clear knee.
