from __future__ import annotations

import csv
import gzip
import json
import math
from dataclasses import dataclass
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "analysis" / "data_runs"
RESULTS_DIR = REPO_ROOT / "analysis" / "results" / "comparison_read_geometry"
NOTEBOOK_PATH = REPO_ROOT / "analysis" / "comparison_read_geometry.ipynb"

INCOMPLETE_LONG_RUN = "run_19_1000_30"
BASELINE_RUN = "run_08_150_30"
PRIMARY_SCOPE = "main_17_common_582"
ALL_RUN_SCOPE = "all_18_common_387"

CLINVAR_GROUPS = ["Benign", "VUS", "Conflicting", "Other", "Pathogenic"]
GOOD_BROAD_GROUPS = {"Benign", "VUS", "Conflicting", "Other"}
GOOD_CORE_GROUPS = {"Benign", "VUS"}


@dataclass(frozen=True)
class RunInfo:
    run_name: str
    run_label: str
    read_len: int
    step: int
    coverage_proxy: float
    variant_rows: int
    successful_genes: int
    failed_genes: int
    incomplete_genes: int
    missing_genes: int


def _read_json(path: Path) -> dict:
    with path.open() as handle:
        return json.load(handle)


def _run_dirs() -> list[Path]:
    return sorted(path.parent for path in DATA_ROOT.glob("run_*/run_manifest.json"))


def load_run_info(run_dir: Path) -> RunInfo:
    manifest = _read_json(run_dir / "run_manifest.json")
    params = _read_json(run_dir / "run_params.json")
    read_generation = params["read_generation"]
    read_len = int(read_generation["read_len"])
    step = int(read_generation["step"])
    counts = manifest["counts"]
    return RunInfo(
        run_name=run_dir.name,
        run_label=str(manifest.get("relative_run_path") or run_dir.name),
        read_len=read_len,
        step=step,
        coverage_proxy=read_len / step,
        variant_rows=int(counts["variant_rows"]),
        successful_genes=int(counts["successful_genes"]),
        failed_genes=int(counts["failed_genes"]),
        incomplete_genes=int(counts["incomplete_genes"]),
        missing_genes=int(counts["missing_genes"]),
    )


def load_gene_rows(run_dir: Path, info: RunInfo) -> list[dict]:
    rows: list[dict] = []
    with gzip.open(run_dir / "genes.csv.gz", "rt", newline="") as handle:
        for row in csv.DictReader(handle):
            rows.append(
                {
                    "run_name": info.run_name,
                    "run_label": info.run_label,
                    "read_len": info.read_len,
                    "step": info.step,
                    "coverage_proxy": info.coverage_proxy,
                    "gene_id": str(row["gene_id"]),
                    "status": row["status"],
                    "total_variants": int(row["total_variants"] or 0),
                    "clinvar_intersected": int(row["clinvar_intersected"] or 0),
                    "clinvar_record_count": int(row["clinvar_record_count"] or 0),
                    "gnomad_intersected": int(row["gnomad_intersected"] or 0),
                    "clinvar_gnomad_both": int(row["clinvar_gnomad_both"] or 0),
                    "mc_annotated_variants": int(row["mc_annotated_variants"] or 0),
                    "csq_annotated_variants": int(row["csq_annotated_variants"] or 0),
                    "error_type": row.get("error_type", ""),
                    "error_message": row.get("error_message", ""),
                }
            )
    return rows


def classify_clnsig(series: pd.Series) -> pd.Series:
    lower = series.fillna("").astype(str).str.lower()
    pathogenic = (
        lower.str.contains("pathogenic", regex=False)
        & ~lower.str.contains("benign", regex=False)
        & ~lower.str.contains("conflicting", regex=False)
    )
    benign = (
        lower.str.contains("benign", regex=False)
        & ~lower.str.contains("pathogenic", regex=False)
        & ~lower.str.contains("conflicting", regex=False)
    )
    vus = lower.str.contains("uncertain", regex=False)
    conflicting = lower.str.contains("conflicting", regex=False) | lower.str.contains(
        "not_provided", regex=False
    )

    grouped = pd.Series("Other", index=series.index, dtype="object")
    grouped.loc[pathogenic] = "Pathogenic"
    grouped.loc[benign] = "Benign"
    grouped.loc[vus] = "VUS"
    grouped.loc[conflicting] = "Conflicting"
    return grouped


def success_gene_sets(gene_df: pd.DataFrame) -> tuple[set[str], set[str]]:
    by_run = {
        run_name: set(rows.loc[rows["status"].eq("success"), "gene_id"].astype(str))
        for run_name, rows in gene_df.groupby("run_name")
    }
    common_all = set.intersection(*by_run.values())
    common_without_long = set.intersection(
        *(genes for run_name, genes in by_run.items() if run_name != INCOMPLETE_LONG_RUN)
    )
    return common_without_long, common_all


def scope_for_run(run_name: str, common_without_long: set[str], common_all: set[str]) -> dict[str, set[str]]:
    scopes = {ALL_RUN_SCOPE: common_all}
    if run_name != INCOMPLETE_LONG_RUN:
        scopes[PRIMARY_SCOPE] = common_without_long
    return scopes


def scan_clinvar_rows(
    run_dirs: list[Path],
    run_info: dict[str, RunInfo],
    common_without_long: set[str],
    common_all: set[str],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    group_counter: dict[tuple[str, str, str], int] = {}
    gene_group_counter: dict[tuple[str, str, str, str], int] = {}
    unique_sets: dict[tuple[str, str, str], set[str]] = {}

    for run_dir in run_dirs:
        info = run_info[run_dir.name]
        scopes = scope_for_run(info.run_name, common_without_long, common_all)
        needed_genes = set.union(*scopes.values())
        usecols = ["gene_id", "variant_key", "clnsig"]

        for chunk in pd.read_csv(
            run_dir / "variants.csv.gz",
            usecols=usecols,
            dtype={"gene_id": "string", "variant_key": "string", "clnsig": "string"},
            chunksize=1_000_000,
        ):
            chunk["clnsig"] = chunk["clnsig"].fillna("")
            chunk = chunk[chunk["clnsig"].ne("")]
            if chunk.empty:
                continue
            chunk = chunk[chunk["gene_id"].isin(needed_genes)]
            if chunk.empty:
                continue
            chunk["clnsig_group"] = classify_clnsig(chunk["clnsig"])
            chunk["gene_variant_key"] = chunk["gene_id"].astype(str) + "|" + chunk["variant_key"].astype(str)

            for scope_name, genes in scopes.items():
                scoped = chunk[chunk["gene_id"].isin(genes)]
                if scoped.empty:
                    continue

                counts = scoped["clnsig_group"].value_counts()
                for group, count in counts.items():
                    group_counter[(scope_name, info.run_name, str(group))] = (
                        group_counter.get((scope_name, info.run_name, str(group)), 0) + int(count)
                    )

                per_gene = scoped.groupby(["gene_id", "clnsig_group"], observed=True).size()
                for (gene_id, group), count in per_gene.items():
                    key = (scope_name, info.run_name, str(gene_id), str(group))
                    gene_group_counter[key] = gene_group_counter.get(key, 0) + int(count)

                for label, groups in [
                    ("clinvar_any", set(CLINVAR_GROUPS)),
                    ("good_broad", GOOD_BROAD_GROUPS),
                    ("good_core", GOOD_CORE_GROUPS),
                    ("pathogenic", {"Pathogenic"}),
                ]:
                    selected = scoped[scoped["clnsig_group"].isin(groups)]
                    if selected.empty:
                        continue
                    for key_type, col in [("variant", "variant_key"), ("gene_variant", "gene_variant_key")]:
                        key = (scope_name, info.run_name, f"{label}_{key_type}")
                        unique_sets.setdefault(key, set()).update(selected[col].dropna().astype(str).tolist())

        print(f"scanned {info.run_name}", flush=True)

    group_rows = [
        {"scope": scope, "run_name": run_name, "clnsig_group": group, "count": count}
        for (scope, run_name, group), count in group_counter.items()
    ]
    gene_group_rows = [
        {
            "scope": scope,
            "run_name": run_name,
            "gene_id": gene_id,
            "clnsig_group": group,
            "count": count,
        }
        for (scope, run_name, gene_id, group), count in gene_group_counter.items()
    ]
    unique_rows = [
        {"scope": scope, "run_name": run_name, "metric": metric, "count": len(values)}
        for (scope, run_name, metric), values in unique_sets.items()
    ]

    unique_overlap_rows = []
    for scope in [PRIMARY_SCOPE, ALL_RUN_SCOPE]:
        baseline_key = (scope, BASELINE_RUN, "good_broad_gene_variant")
        if baseline_key not in unique_sets:
            continue
        baseline = unique_sets[baseline_key]
        for run_name in sorted({key[1] for key in unique_sets if key[0] == scope}):
            run_set = unique_sets.get((scope, run_name, "good_broad_gene_variant"), set())
            unique_overlap_rows.append(
                {
                    "scope": scope,
                    "run_name": run_name,
                    "baseline_run": BASELINE_RUN,
                    "good_broad_gene_variant": len(run_set),
                    "new_vs_baseline": len(run_set - baseline),
                    "lost_vs_baseline": len(baseline - run_set),
                    "shared_with_baseline": len(run_set & baseline),
                    "jaccard_vs_baseline": len(run_set & baseline) / len(run_set | baseline)
                    if (run_set | baseline)
                    else math.nan,
                }
            )

    return (
        pd.DataFrame(group_rows),
        pd.DataFrame(gene_group_rows),
        pd.DataFrame(unique_rows + unique_overlap_rows),
    )


def build_scope_metrics(
    gene_df: pd.DataFrame,
    group_df: pd.DataFrame,
    unique_df: pd.DataFrame,
    run_info: dict[str, RunInfo],
    common_without_long: set[str],
    common_all: set[str],
) -> pd.DataFrame:
    rows = []
    for run_name, info in run_info.items():
        for scope, genes in scope_for_run(run_name, common_without_long, common_all).items():
            scoped_genes = gene_df[
                gene_df["run_name"].eq(run_name)
                & gene_df["gene_id"].isin(genes)
                & gene_df["status"].eq("success")
            ]
            class_counts = {group: 0 for group in CLINVAR_GROUPS}
            scoped_groups = group_df[group_df["scope"].eq(scope) & group_df["run_name"].eq(run_name)]
            for _, row in scoped_groups.iterrows():
                class_counts[str(row["clnsig_group"])] = int(row["count"])

            scoped_unique = unique_df[
                unique_df["scope"].eq(scope)
                & unique_df["run_name"].eq(run_name)
                & unique_df["metric"].notna()
            ]
            unique_counts = {
                str(row["metric"]): int(row["count"])
                for _, row in scoped_unique.iterrows()
                if str(row["metric"]).endswith("_variant") or str(row["metric"]).endswith("_gene_variant")
            }

            total_variants = int(scoped_genes["total_variants"].sum())
            clinvar_any = sum(class_counts.values())
            good_broad = sum(class_counts[group] for group in GOOD_BROAD_GROUPS)
            good_core = sum(class_counts[group] for group in GOOD_CORE_GROUPS)
            pathogenic = class_counts["Pathogenic"]
            rows.append(
                {
                    "scope": scope,
                    "run_name": run_name,
                    "read_len": info.read_len,
                    "step": info.step,
                    "geometry": f"{info.read_len}/{info.step}",
                    "coverage_proxy": info.coverage_proxy,
                    "successful_genes_in_scope": len(scoped_genes),
                    "total_variants": total_variants,
                    "clinvar_any": clinvar_any,
                    "good_broad": good_broad,
                    "good_core": good_core,
                    "pathogenic": pathogenic,
                    "benign": class_counts["Benign"],
                    "vus": class_counts["VUS"],
                    "conflicting": class_counts["Conflicting"],
                    "other_clinvar": class_counts["Other"],
                    "good_broad_share_pct": 100 * good_broad / total_variants if total_variants else math.nan,
                    "good_core_share_pct": 100 * good_core / total_variants if total_variants else math.nan,
                    "clinvar_share_pct": 100 * clinvar_any / total_variants if total_variants else math.nan,
                    "variants_per_good_broad": total_variants / good_broad if good_broad else math.nan,
                    "nonclinvar_per_good_broad": (total_variants - clinvar_any) / good_broad if good_broad else math.nan,
                    "unique_good_broad_variants": unique_counts.get("good_broad_variant", 0),
                    "unique_good_core_variants": unique_counts.get("good_core_variant", 0),
                    "unique_good_broad_gene_variants": unique_counts.get("good_broad_gene_variant", 0),
                    "unique_good_core_gene_variants": unique_counts.get("good_core_gene_variant", 0),
                    "unique_pathogenic_gene_variants": unique_counts.get("pathogenic_gene_variant", 0),
                }
            )
    return pd.DataFrame(rows)


def add_decision_metrics(metrics: pd.DataFrame) -> pd.DataFrame:
    out = metrics.copy()
    out["good_broad_norm"] = out.groupby("scope")["good_broad"].transform(_minmax)
    out["total_variants_norm"] = out.groupby("scope")["total_variants"].transform(_minmax)
    out["coverage_norm"] = out.groupby("scope")["coverage_proxy"].transform(_minmax)
    out["score_lambda_025"] = out["good_broad_norm"] - 0.25 * out["total_variants_norm"]
    out["score_lambda_050"] = out["good_broad_norm"] - 0.50 * out["total_variants_norm"]
    out["score_lambda_075"] = out["good_broad_norm"] - 0.75 * out["total_variants_norm"]

    knee_rows = []
    for scope, scoped in out.groupby("scope"):
        scoped = scoped.sort_values("total_variants").reset_index(drop=True)
        x = scoped["total_variants"].astype(float)
        y = scoped["good_broad"].astype(float)
        x_norm = (x - x.min()) / (x.max() - x.min()) if x.max() != x.min() else x * 0
        y_norm = (y - y.min()) / (y.max() - y.min()) if y.max() != y.min() else y * 0
        # Distance from the diagonal joining the least/most costly points.
        distance = (y_norm - x_norm) / math.sqrt(2)
        for run_name, dist in zip(scoped["run_name"], distance):
            knee_rows.append({"scope": scope, "run_name": run_name, "knee_distance": float(dist)})
    knee_df = pd.DataFrame(knee_rows)
    out = out.merge(knee_df, on=["scope", "run_name"], how="left")
    out["is_pareto_good_vs_variants"] = False
    for scope, scoped in out.groupby("scope"):
        idxs = []
        for idx, row in scoped.iterrows():
            dominated = scoped[
                (scoped["good_broad"] >= row["good_broad"])
                & (scoped["total_variants"] <= row["total_variants"])
                & (
                    (scoped["good_broad"] > row["good_broad"])
                    | (scoped["total_variants"] < row["total_variants"])
                )
            ]
            if dominated.empty:
                idxs.append(idx)
        out.loc[idxs, "is_pareto_good_vs_variants"] = True
    return out


def _minmax(series: pd.Series) -> pd.Series:
    if series.max() == series.min():
        return pd.Series(0.0, index=series.index)
    return (series - series.min()) / (series.max() - series.min())


def build_increment_table(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scope, scoped in metrics.groupby("scope"):
        scoped = scoped.sort_values("total_variants")
        prev = None
        for _, row in scoped.iterrows():
            if prev is None:
                delta_good = math.nan
                delta_variants = math.nan
                marginal = math.nan
            else:
                delta_good = int(row["good_broad"] - prev["good_broad"])
                delta_variants = int(row["total_variants"] - prev["total_variants"])
                marginal = 1_000_000 * delta_good / delta_variants if delta_variants else math.nan
            rows.append(
                {
                    "scope": scope,
                    "run_name": row["run_name"],
                    "geometry": row["geometry"],
                    "total_variants": int(row["total_variants"]),
                    "good_broad": int(row["good_broad"]),
                    "delta_good_vs_previous_cost": delta_good,
                    "delta_variants_vs_previous_cost": delta_variants,
                    "marginal_good_per_million_variants": marginal,
                }
            )
            prev = row
    return pd.DataFrame(rows)


def write_plots(metrics: pd.DataFrame, overlap: pd.DataFrame) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    main = metrics[metrics["scope"].eq(PRIMARY_SCOPE)].copy()
    main = main.sort_values("coverage_proxy")

    plt.figure(figsize=(10, 6))
    plt.scatter(main["total_variants"] / 1e6, main["good_broad"] / 1000, s=70)
    for _, row in main.iterrows():
        plt.annotate(row["geometry"], (row["total_variants"] / 1e6, row["good_broad"] / 1000), fontsize=8)
    plt.xlabel("Total variants, M")
    plt.ylabel("Good ClinVar non-pathogenic/VUS rows, K")
    plt.title("Good ClinVar capture vs downstream variant burden")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "good_vs_total_variants.png", dpi=180)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.plot(main["coverage_proxy"], main["good_broad"], marker="o", label="Good broad")
    plt.plot(main["coverage_proxy"], main["good_core"], marker="o", label="Benign + VUS")
    for _, row in main.iterrows():
        plt.annotate(row["geometry"], (row["coverage_proxy"], row["good_broad"]), fontsize=8)
    plt.xlabel("Coverage proxy: read_len / step")
    plt.ylabel("ClinVar rows")
    plt.title("Saturation by coverage proxy")
    plt.grid(alpha=0.25)
    plt.legend()
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "saturation_by_coverage_proxy.png", dpi=180)
    plt.close()

    plt.figure(figsize=(10, 6))
    plt.scatter(main["good_broad"], main["good_broad_share_pct"], s=70)
    for _, row in main.iterrows():
        plt.annotate(row["geometry"], (row["good_broad"], row["good_broad_share_pct"]), fontsize=8)
    plt.xlabel("Good ClinVar non-pathogenic/VUS rows")
    plt.ylabel("Good share of total variants, %")
    plt.title("Absolute gain comes with falling ClinVar density")
    plt.grid(alpha=0.25)
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / "good_count_vs_good_share.png", dpi=180)
    plt.close()

    overlap_main = overlap[overlap["scope"].eq(PRIMARY_SCOPE)].copy()
    if not overlap_main.empty:
        overlap_main = overlap_main.merge(main[["run_name", "geometry"]], on="run_name", how="left")
        overlap_main = overlap_main.sort_values("new_vs_baseline", ascending=False)
        plt.figure(figsize=(11, 6))
        x = range(len(overlap_main))
        plt.bar(x, overlap_main["new_vs_baseline"], label="New vs 150/30")
        plt.bar(x, -overlap_main["lost_vs_baseline"], label="Lost vs 150/30")
        plt.xticks(list(x), overlap_main["geometry"], rotation=45, ha="right")
        plt.ylabel("Good gene-variant rows")
        plt.title("Good ClinVar set movement relative to 150/30")
        plt.grid(axis="y", alpha=0.25)
        plt.legend()
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / "baseline_overlap_150_30.png", dpi=180)
        plt.close()


def write_notebook() -> None:
    def md(source: str) -> dict:
        return {"cell_type": "markdown", "metadata": {}, "source": source.strip().splitlines(True)}

    def code(source: str) -> dict:
        return {
            "cell_type": "code",
            "execution_count": None,
            "metadata": {},
            "outputs": [],
            "source": source.strip().splitlines(True),
        }

    notebook = {
        "cells": [
            md(
                """
                # Comparison Read Geometry

                Цель: выбрать read geometry, которая максимизирует абсолютный охват ClinVar вариантов,
                не являющихся явно pathogenic/likely pathogenic, и при этом увидеть насыщение по цене
                общего числа variant calls.

                Основной анализ использует 17 завершенных run'ов на общем наборе 582 успешных генов.
                `1000/30` анализируется отдельно на общем пересечении 387 генов, потому что экспорт этого
                run'а неполный: 387 success, 21 incomplete, 178 missing.
                """
            ),
            code(
                """
                from pathlib import Path
                import pandas as pd

                RESULTS = Path('results/comparison_read_geometry')
                if not RESULTS.exists():
                    RESULTS = Path('analysis/results/comparison_read_geometry')

                metrics = pd.read_csv(RESULTS / 'scope_metrics.csv')
                increments = pd.read_csv(RESULTS / 'incremental_cost_curve.csv')
                overlap = pd.read_csv(RESULTS / 'baseline_overlap.csv')
                per_gene = pd.read_csv(RESULTS / 'per_gene_metrics.csv')

                main = metrics[metrics.scope.eq('main_17_common_582')].sort_values('good_broad', ascending=False)
                main[['run_name','geometry','coverage_proxy','total_variants','good_broad','good_core','pathogenic',
                      'good_broad_share_pct','variants_per_good_broad','knee_distance',
                      'score_lambda_050','is_pareto_good_vs_variants']].head(20)
                """
            ),
            md(
                """
                ## Definition of "good"

                Primary metric `good_broad` = ClinVar rows grouped as Benign, VUS, Conflicting, or Other.
                Rows grouped as Pathogenic are excluded. The conservative sensitivity metric `good_core`
                keeps only Benign + VUS.

                This keeps the selection aligned with the project goal: catch as many non-pathogenic/VUS
                ClinVar-supported variants as possible, while postponing stricter filtering to later stages.
                """
            ),
            md(
                """
                ## Decision Method

                I use three complementary criteria rather than a single raw count:

                1. **Absolute capture:** maximize `good_broad`.
                2. **Saturation / burden:** track `good_broad_share_pct` and `variants_per_good_broad`;
                   if good variants grow slowly while total calls explode, the configuration is past the shoulder.
                3. **Dominance:** reject geometries that produce fewer good variants and more total variants than
                   another geometry. This removes `150/20`, `150/25`, and `1000/60` from the main candidate set.

                A strict max-capture choice is `500/30`. A pragmatic max-capture choice is `300/30`: it captures
                130,771 good rows, which is 98.4% of the `500/30` maximum, while avoiding the longest/highest-burden
                setting. The statistical elbow is around `225/45` / `200/40`, where the good count is already about
                92% of the observed maximum with much lower total variant burden.

                | Geometry | Good broad | Total variants | Good share | Role |
                |---:|---:|---:|---:|---|
                | 500/30 | 132,847 | 23.01M | 0.577% | observed max |
                | 300/30 | 130,771 | 21.67M | 0.603% | pragmatic max-capture choice |
                | 225/45 | 122,583 | 14.28M | 0.858% | elbow / balanced choice |
                | 200/40 | 121,577 | 14.02M | 0.867% | similar elbow, slightly lower capture |
                | 150/30 | 118,365 | 13.22M | 0.895% | current baseline |
                """
            ),
            code(
                """
                main[['geometry','benign','vus','conflicting','other_clinvar','pathogenic','good_broad','good_core']]
                """
            ),
            md("## Saturation Figures\n\n![Good vs total variants](results/comparison_read_geometry/good_vs_total_variants.png)"),
            md("![Saturation by coverage proxy](results/comparison_read_geometry/saturation_by_coverage_proxy.png)"),
            md("![Good count vs good share](results/comparison_read_geometry/good_count_vs_good_share.png)"),
            md("![Baseline overlap](results/comparison_read_geometry/baseline_overlap_150_30.png)"),
            md(
                """
                ## Incremental Cost Curve

                The table is ordered by observed total variant burden, not by run number. The marginal column
                estimates how many additional good ClinVar rows are gained per additional million total variant
                calls when moving to the next more expensive observed point.
                """
            ),
            code(
                """
                increments[increments.scope.eq('main_17_common_582')].sort_values('total_variants')
                """
            ),
            md(
                """
                ## 1000/30 Sensitivity

                Because `1000/30` is incomplete, compare it only on the 387 genes that succeeded in all 18 runs.
                On that restricted set it still does not beat `500/30`, and it is effectively tied with `300/30`
                / `1000/60` while being operationally much worse. I would not spend more cluster time on completing
                `1000/30` unless the project explicitly needs the extreme endpoint.
                """
            ),
            code(
                """
                metrics[metrics.scope.eq('all_18_common_387')].sort_values('good_broad', ascending=False)[
                    ['run_name','geometry','successful_genes_in_scope','total_variants','good_broad','good_core',
                     'good_broad_share_pct','variants_per_good_broad','knee_distance']
                ]
                """
            ),
            md(
                """
                ## Recommended Next Runs

                If more runs are needed, keep them few and use them to separate read-length effects from coverage
                effects without going back to `1000/30`:

                - `400/40`: same coverage proxy as `300/30` (`read_len / step = 10`) with longer reads. Tests whether
                  the 300 -> 500 gain is mainly read length or just extra coverage.
                - `500/50`: same 500 bp read length as the current max, but coverage proxy 10 instead of 16.7.
                  Tests whether `500/30` can be made much cheaper without losing most good variants.
                - Optional only if the first two are ambiguous: `300/40`, a cheaper 300 bp setting that tests whether
                  the chosen default can be relaxed from step 30.

                I would stop running denser 150 bp geometries: `150/20` is dominated by `200/30`, and the 150 bp
                step sweep has already shown diminishing returns.
                """
            ),
            md(
                """
                ## Re-run Heavy Scan

                The notebook reads precomputed tables. To refresh them after new exports:

                ```bash
                MPLCONFIGDIR=/tmp/matplotlib-cache .venv/bin/python3 analysis/comparison_read_geometry_analysis.py
                ```
                """
            ),
        ],
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
        },
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    NOTEBOOK_PATH.write_text(json.dumps(notebook, indent=2, ensure_ascii=False) + "\n")


def main() -> int:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    run_dirs = _run_dirs()
    infos = {run_dir.name: load_run_info(run_dir) for run_dir in run_dirs}

    gene_df = pd.DataFrame(
        [row for run_dir in run_dirs for row in load_gene_rows(run_dir, infos[run_dir.name])]
    )
    common_without_long, common_all = success_gene_sets(gene_df)
    gene_df.to_csv(RESULTS_DIR / "per_gene_metrics.csv", index=False)

    run_overview = pd.DataFrame([info.__dict__ for info in infos.values()]).sort_values(
        ["read_len", "step", "run_name"]
    )
    run_overview["primary_common_genes"] = len(common_without_long)
    run_overview["all_run_common_genes"] = len(common_all)
    run_overview.to_csv(RESULTS_DIR / "run_overview.csv", index=False)

    group_df, gene_group_df, unique_df = scan_clinvar_rows(
        run_dirs=run_dirs,
        run_info=infos,
        common_without_long=common_without_long,
        common_all=common_all,
    )
    group_df.to_csv(RESULTS_DIR / "clinvar_group_counts.csv", index=False)
    gene_group_df.to_csv(RESULTS_DIR / "per_gene_clinvar_group_counts.csv", index=False)

    baseline_overlap = unique_df[unique_df["baseline_run"].notna()].copy()
    baseline_overlap.to_csv(RESULTS_DIR / "baseline_overlap.csv", index=False)
    unique_df[unique_df["baseline_run"].isna()].to_csv(RESULTS_DIR / "unique_clinvar_counts.csv", index=False)

    metrics = build_scope_metrics(gene_df, group_df, unique_df, infos, common_without_long, common_all)
    metrics = add_decision_metrics(metrics)
    metrics = metrics.sort_values(["scope", "good_broad"], ascending=[True, False])
    metrics.to_csv(RESULTS_DIR / "scope_metrics.csv", index=False)

    increments = build_increment_table(metrics)
    increments.to_csv(RESULTS_DIR / "incremental_cost_curve.csv", index=False)

    write_plots(metrics, baseline_overlap)
    write_notebook()

    print(f"Wrote results to {RESULTS_DIR}")
    print(f"Wrote notebook to {NOTEBOOK_PATH}")
    print(f"Primary common genes: {len(common_without_long)}")
    print(f"All-run common genes: {len(common_all)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
