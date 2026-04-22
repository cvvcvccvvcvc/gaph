from __future__ import annotations

import argparse
import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from run_snapshot import discover_exported_run_dirs, load_exported_run_data


IMPACT_RANK = {"MODIFIER": 0, "LOW": 1, "MODERATE": 2, "HIGH": 3, "OTHER": -1}

TERM_TO_IMPACT = {
    "frameshift_variant": "HIGH",
    "stop_gained": "HIGH",
    "stop_lost": "HIGH",
    "start_lost": "HIGH",
    "splice_acceptor_variant": "HIGH",
    "splice_donor_variant": "HIGH",
    "missense_variant": "MODERATE",
    "inframe_insertion": "MODERATE",
    "inframe_deletion": "MODERATE",
    "protein_altering_variant": "MODERATE",
    "splice_region_variant": "MODERATE",
    "synonymous_variant": "LOW",
    "intron_variant": "MODIFIER",
    "5_prime_utr_variant": "MODIFIER",
    "3_prime_utr_variant": "MODIFIER",
    "upstream_gene_variant": "MODIFIER",
    "downstream_gene_variant": "MODIFIER",
    "intergenic_variant": "MODIFIER",
    "regulatory_region_variant": "MODIFIER",
    "non_coding_transcript_variant": "MODIFIER",
    "coding_sequence_variant": "MODIFIER",
    "mature_mirna_variant": "MODIFIER",
}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare compact run snapshots exported into analysis/data_runs."
    )
    parser.add_argument(
        "exports_root",
        type=Path,
        help="Root directory that contains exported run snapshots.",
    )
    parser.add_argument(
        "--run",
        action="append",
        default=[],
        help="Run directory name or relative path under exports_root. Can be passed multiple times.",
    )
    parser.add_argument(
        "--gene",
        action="append",
        default=[],
        help="Gene ID to include. Can be passed multiple times.",
    )
    parser.add_argument(
        "--genes-file",
        type=Path,
        help="Optional file with one gene ID per line.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("analysis/results/csv_compare"),
        help="Directory for comparison CSV outputs.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    genes_filter = _load_gene_filter(args.gene, args.genes_file)
    run_dirs = _select_run_dirs(args.exports_root, args.run)
    if not run_dirs:
        print(f"No exported runs found under {args.exports_root}")
        return 1

    run_data = {}
    include_variants = genes_filter is not None
    for export_dir in run_dirs:
        run_label = export_dir.relative_to(args.exports_root).as_posix()
        if run_label == ".":
            run_label = export_dir.name
        run_data[run_label] = load_exported_run_data(
            export_dir,
            genes_filter=genes_filter,
            run_label=run_label,
            include_variants=include_variants,
        )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    outputs = build_outputs(run_data, args.output_dir)

    print(f"Loaded {len(run_data)} runs from {args.exports_root}")
    print("\nOverall statistics:")
    print(outputs["df_stats"].to_string(index=False))
    print(f"\nSaved comparison tables to {args.output_dir}")
    return 0


def build_outputs(run_data: dict[str, dict], output_dir: Path) -> dict[str, pd.DataFrame]:
    output_dir.mkdir(parents=True, exist_ok=True)
    run_ids = list(run_data.keys())

    df_params = _build_params_comparison(run_data)
    df_stats = _build_overall_stats(run_data)
    df_clnsig = _build_clnsig_distribution(run_data)
    df_grouped = _build_grouped_clnsig(run_data)
    df_benign_stars = _build_benign_stars(run_data)
    df_conflicting_dist, conflicting_gene_tables = _build_conflicting_tables(run_data)
    df_af_stats, df_af_categories = _build_af_tables(run_data)
    consequence_outputs = _build_consequence_outputs(run_data)

    df_params.to_csv(output_dir / "params_comparison.csv", index=False)
    df_stats.to_csv(output_dir / "overall_statistics.csv", index=False)
    df_clnsig.to_csv(output_dir / "clnsig_distribution.csv", index=False)
    df_grouped.to_csv(output_dir / "grouped_clnsig.csv", index=False)
    df_af_stats.to_csv(output_dir / "af_statistics.csv", index=False)
    df_af_categories.to_csv(output_dir / "af_categories.csv", index=False)

    if not df_benign_stars.empty:
        df_benign_stars.to_csv(output_dir / "benign_stars_by_run.csv", index=False)
    if not df_conflicting_dist.empty:
        df_conflicting_dist.to_csv(output_dir / "conflicting_clnsigconf_distribution.csv", index=False)
    for run_id, df_gene_conf in conflicting_gene_tables.items():
        if not df_gene_conf.empty:
            safe_run_id = run_id.replace("/", "__")
            df_gene_conf.to_csv(output_dir / f"conflicting_clnsigconf_by_gene_{safe_run_id}.csv", index=False)

    for name, df in consequence_outputs.items():
        if isinstance(df, pd.DataFrame) and not df.empty:
            df.to_csv(output_dir / f"{name}.csv", index=False)

    return {
        "df_params": df_params,
        "df_stats": df_stats,
        "df_clnsig": df_clnsig,
        "df_grouped": df_grouped,
        "df_benign_stars": df_benign_stars,
        "df_conflicting_dist": df_conflicting_dist,
        "df_af_stats": df_af_stats,
        "df_af_categories": df_af_categories,
        **consequence_outputs,
    }


def _load_gene_filter(genes: list[str], genes_file: Path | None) -> set[str] | None:
    gene_ids = {str(gene).strip() for gene in genes if str(gene).strip()}
    if genes_file:
        with genes_file.open() as handle:
            for line in handle:
                line = line.strip()
                if line:
                    gene_ids.add(line)
    return gene_ids or None


def _select_run_dirs(exports_root: Path, selected_runs: list[str]) -> list[Path]:
    run_dirs = discover_exported_run_dirs(exports_root)
    if not selected_runs:
        return run_dirs

    selected = set(selected_runs)
    return [
        run_dir
        for run_dir in run_dirs
        if run_dir.name in selected or run_dir.relative_to(exports_root).as_posix() in selected
    ]


def _build_params_comparison(run_data: dict[str, dict]) -> pd.DataFrame:
    params_rows = []
    all_param_keys = set()

    for data in run_data.values():
        all_param_keys.update(data["params"].keys())

    ignored = {
        "gene_ids",
        "ortholog_source",
        "hitlist_size",
        "cache",
        "blast_expect",
        "config_path",
        "run_dir",
        "started_at",
    }

    for key in sorted(all_param_keys):
        if key in ignored:
            continue
        row = {"Parameter": key}
        for run_id, data in run_data.items():
            value = data["params"].get(key, "N/A")
            row[run_id] = str(value) if isinstance(value, list) else value
        params_rows.append(row)

    return pd.DataFrame(params_rows)


def _build_overall_stats(run_data: dict[str, dict]) -> pd.DataFrame:
    stats_rows = []
    for run_id, data in run_data.items():
        gene_results = data["gene_results"]
        total_genes = len(gene_results)
        total_variants = sum(g["total_variants"] for g in gene_results)
        total_clinvar = sum(g["clinvar_intersected"] for g in gene_results)
        total_gnomad = sum(g["gnomad_intersected"] for g in gene_results)
        total_both = sum(g["clinvar_gnomad_both"] for g in gene_results)

        clinvar_pct = (total_clinvar / total_variants * 100) if total_variants else 0
        gnomad_pct = (total_gnomad / total_variants * 100) if total_variants else 0
        both_pct = (total_both / total_variants * 100) if total_variants else 0

        stats_rows.append(
            {
                "Run": run_id,
                "Genes": total_genes,
                "Total Variants": f"{total_variants:,}",
                "ClinVar": f"{total_clinvar:,} ({clinvar_pct:.2f}%)",
                "gnomAD": f"{total_gnomad:,} ({gnomad_pct:.2f}%)",
                "ClinVar+gnomAD": f"{total_both:,} ({both_pct:.2f}%)",
                "_total_variants_numeric": total_variants,
            }
        )

    df_stats = pd.DataFrame(stats_rows)
    if not df_stats.empty:
        df_stats = df_stats.sort_values("_total_variants_numeric", ascending=False).drop(columns=["_total_variants_numeric"])
    return df_stats


def _build_clnsig_distribution(run_data: dict[str, dict]) -> pd.DataFrame:
    clnsig_by_run = {}
    for run_id, data in run_data.items():
        summary = data.get("analysis_summary")
        if summary:
            clnsig_by_run[run_id] = _summary_count_map(summary.get("clnsig_counts", []), "CLNSIG")
        else:
            counts = defaultdict(int)
            for clnsig in data["all_clnsig"]:
                for sig in split_clnsig_labels(clnsig):
                    counts[sig] += 1
            clnsig_by_run[run_id] = counts

    all_clnsig_types = set()
    for counts in clnsig_by_run.values():
        all_clnsig_types.update(counts.keys())

    rows = []
    run_ids = list(run_data.keys())
    for clnsig_type in sorted(all_clnsig_types):
        row = {"CLNSIG": clnsig_type}
        for run_id in run_ids:
            row[run_id] = clnsig_by_run[run_id].get(clnsig_type, 0)
        rows.append(row)

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    if len(run_ids) == 2:
        df["Δ"] = df[run_ids[1]] - df[run_ids[0]]
    df["_total"] = df[run_ids].sum(axis=1)
    return df.sort_values("_total", ascending=False).drop(columns=["_total"])


def _build_grouped_clnsig(run_data: dict[str, dict]) -> pd.DataFrame:
    grouped_by_run = {}
    for run_id, data in run_data.items():
        summary = data.get("analysis_summary")
        if summary:
            grouped_by_run[run_id] = _summary_count_map(summary.get("grouped_clnsig_counts", []), "Group")
        else:
            counts = defaultdict(int)
            for clnsig in data["all_clnsig"]:
                counts[group_clnsig(clnsig)] += 1
            grouped_by_run[run_id] = counts

    rows = []
    group_order = ["Benign", "Pathogenic", "VUS", "Conflicting", "Other"]
    run_ids = list(run_data.keys())
    for group in group_order:
        row = {"Group": group}
        for run_id in run_ids:
            row[run_id] = grouped_by_run[run_id].get(group, 0)
        if any(row[run_id] > 0 for run_id in run_ids):
            rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty and len(run_ids) == 2:
        df["Δ"] = df[run_ids[1]] - df[run_ids[0]]
    return df


def _build_benign_stars(run_data: dict[str, dict]) -> pd.DataFrame:
    benign_stars_by_run = {}
    for run_id, data in run_data.items():
        summary = data.get("analysis_summary")
        if summary:
            benign_stars_by_run[run_id] = _summary_count_map(summary.get("benign_star_counts", []), "benign_star_bucket")
        else:
            counts = defaultdict(int)
            for rec in data["all_clinvar_records"]:
                labels = split_clnsig_labels(rec.get("clnsig", ""))
                if not labels:
                    continue
                if any(is_benign_label(label) for label in labels):
                    stars = clinvar_stars_from_review_status(rec.get("clnrevstat", ""))
                    bucket = f"benign_{stars}" if stars is not None else "benign_unknown"
                    counts[bucket] += 1
            benign_stars_by_run[run_id] = counts

    run_ids = list(run_data.keys())
    star_order = [f"benign_{i}" for i in range(5)] + ["benign_unknown"]
    rows = []
    for bucket in star_order:
        row = {"benign_star_bucket": bucket}
        for run_id in run_ids:
            row[run_id] = benign_stars_by_run[run_id].get(bucket, 0)
        rows.append(row)

    df = pd.DataFrame(rows)
    if not df.empty and len(run_ids) == 2:
        df["Δ"] = df[run_ids[1]] - df[run_ids[0]]
    return df


def _build_conflicting_tables(run_data: dict[str, dict]) -> tuple[pd.DataFrame, dict[str, pd.DataFrame]]:
    run_ids = list(run_data.keys())
    conflicting_gene_tables = {}
    conflicting_submitter_dist_by_run = {}
    all_conf_labels = set()
    run_gene_maps = {}

    for run_id, data in run_data.items():
        summary = data.get("analysis_summary")
        if summary:
            summary_counts = _summary_count_map(summary.get("conflicting_clnsigconf_counts", []), "CLNSIGCONF")
            conflicting_submitter_dist_by_run[run_id] = defaultdict(int, summary_counts)
            all_conf_labels.update(summary_counts.keys())

            gene_rows = summary.get("conflicting_clnsigconf_by_gene", [])
            if gene_rows:
                conflicting_gene_tables[run_id] = (
                    pd.DataFrame(gene_rows)
                    .sort_values("total_submitters", ascending=False)
                    .reset_index(drop=True)
                )
            else:
                conflicting_gene_tables[run_id] = pd.DataFrame()
            continue

        gene_map = {}
        submitter_dist = defaultdict(int)

        for gene in data["gene_results"]:
            gene_id = str(gene["gene_id"])
            per_gene = defaultdict(int)

            for rec in gene.get("clinvar_records", []):
                labels = split_clnsig_labels(rec.get("clnsig", ""))
                if "Conflicting_classifications_of_pathogenicity" not in labels:
                    continue
                parsed = parse_clnsigconf_submitter_counts(rec.get("clnsigconf", ""))
                for label, count in parsed.items():
                    per_gene[label] += int(count)
                    submitter_dist[label] += int(count)
                    all_conf_labels.add(label)

            if per_gene:
                gene_map[gene_id] = per_gene

        run_gene_maps[run_id] = gene_map
        conflicting_submitter_dist_by_run[run_id] = submitter_dist

    conf_labels_sorted = sorted(all_conf_labels, key=lambda item: item.lower())
    for run_id in run_ids:
        gene_map = run_gene_maps.get(run_id, {})
        if not gene_map:
            conflicting_gene_tables[run_id] = pd.DataFrame()
            continue

        rows = []
        for gene_id in sorted(gene_map.keys(), key=_gene_sort_key):
            row = {"gene_id": gene_id}
            per_gene = gene_map[gene_id]
            for label in conf_labels_sorted:
                row[label] = int(per_gene.get(label, 0))
            row["total_submitters"] = int(sum(per_gene.values()))
            rows.append(row)
        conflicting_gene_tables[run_id] = (
            pd.DataFrame(rows)
            .sort_values("total_submitters", ascending=False)
            .reset_index(drop=True)
        )

    if not conf_labels_sorted:
        return pd.DataFrame(), conflicting_gene_tables

    rows = []
    for label in conf_labels_sorted:
        row = {"CLNSIGCONF": label}
        for run_id in run_ids:
            row[run_id] = int(conflicting_submitter_dist_by_run[run_id].get(label, 0))
        rows.append(row)

    df = pd.DataFrame(rows)
    if len(run_ids) == 2:
        df["Δ"] = df[run_ids[1]] - df[run_ids[0]]
    df["_total"] = df[run_ids].sum(axis=1)
    df = df.sort_values("_total", ascending=False).drop(columns=["_total"])
    return df, conflicting_gene_tables


def _build_af_tables(run_data: dict[str, dict]) -> tuple[pd.DataFrame, pd.DataFrame]:
    af_stats_rows = []
    af_category_rows = []

    for run_id, data in run_data.items():
        summary = data.get("analysis_summary")
        if summary:
            af_stats = summary.get("af_stats", {})
            count = int(af_stats.get("count") or 0)
            if count:
                rare = int(af_stats.get("rare", 0))
                low_freq = int(af_stats.get("low_freq", 0))
                common = int(af_stats.get("common", 0))
                af_stats_rows.append(
                    {
                        "Run": run_id,
                        "Count": count,
                        "Min": f"{float(af_stats['min']):.6f}",
                        "Max": f"{float(af_stats['max']):.6f}",
                        "Mean": f"{float(af_stats['mean']):.6f}",
                        "Median": f"{float(af_stats['median']):.6f}",
                    }
                )
                af_category_rows.append(
                    {
                        "Run": run_id,
                        "Rare (<1%)": f"{rare} ({rare/count*100:.1f}%)",
                        "Low freq (1-5%)": f"{low_freq} ({low_freq/count*100:.1f}%)",
                        "Common (>=5%)": f"{common} ({common/count*100:.1f}%)",
                    }
                )
            else:
                af_stats_rows.append(
                    {"Run": run_id, "Count": 0, "Min": "N/A", "Max": "N/A", "Mean": "N/A", "Median": "N/A"}
                )
                af_category_rows.append(
                    {"Run": run_id, "Rare (<1%)": "N/A", "Low freq (1-5%)": "N/A", "Common (>=5%)": "N/A"}
                )
            continue

        af_values = data["all_gnomad_af"]
        if af_values:
            af_series = pd.Series(af_values)
            rare = (af_series < 0.01).sum()
            low_freq = ((af_series >= 0.01) & (af_series < 0.05)).sum()
            common = (af_series >= 0.05).sum()

            af_stats_rows.append(
                {
                    "Run": run_id,
                    "Count": len(af_values),
                    "Min": f"{af_series.min():.6f}",
                    "Max": f"{af_series.max():.6f}",
                    "Mean": f"{af_series.mean():.6f}",
                    "Median": f"{af_series.median():.6f}",
                }
            )
            af_category_rows.append(
                {
                    "Run": run_id,
                    "Rare (<1%)": f"{rare} ({rare/len(af_values)*100:.1f}%)",
                    "Low freq (1-5%)": f"{low_freq} ({low_freq/len(af_values)*100:.1f}%)",
                    "Common (>=5%)": f"{common} ({common/len(af_values)*100:.1f}%)",
                }
            )
        else:
            af_stats_rows.append(
                {"Run": run_id, "Count": 0, "Min": "N/A", "Max": "N/A", "Mean": "N/A", "Median": "N/A"}
            )
            af_category_rows.append(
                {"Run": run_id, "Rare (<1%)": "N/A", "Low freq (1-5%)": "N/A", "Common (>=5%)": "N/A"}
            )

    return pd.DataFrame(af_stats_rows), pd.DataFrame(af_category_rows)


def _build_consequence_outputs(run_data: dict[str, dict]) -> dict[str, pd.DataFrame]:
    if all(data.get("analysis_summary") for data in run_data.values()):
        run_impact_rows = []
        run_gene_impact_rows = []
        run_mismatch_rows = []
        run_gene_mismatch_rows = []
        mc_term_rows = []
        csq_term_rows = []

        for run_id, data in run_data.items():
            summary = data["analysis_summary"]
            run_impact_rows.extend([{**row, "run_id": run_id} for row in summary.get("run_impact_rows", [])])
            run_gene_impact_rows.extend([{**row, "run_id": run_id} for row in summary.get("run_gene_impact_rows", [])])
            if summary.get("run_mismatch_summary"):
                run_mismatch_rows.append({"run_id": run_id, **summary["run_mismatch_summary"]})
            run_gene_mismatch_rows.extend([{**row, "run_id": run_id} for row in summary.get("run_gene_mismatch_rows", [])])
            mc_term_rows.extend([{"run_id": run_id, **row} for row in summary.get("mc_terms", [])])
            csq_term_rows.extend([{"run_id": run_id, **row} for row in summary.get("csq_terms", [])])

        return {
            "run_impact_class_shares": pd.DataFrame(run_impact_rows),
            "run_gene_impact_class_shares": pd.DataFrame(run_gene_impact_rows),
            "run_mc_vs_csq_mismatch_summary": pd.DataFrame(run_mismatch_rows),
            "run_gene_mc_vs_csq_mismatch": pd.DataFrame(run_gene_mismatch_rows),
            "run_mc_terms": pd.DataFrame(mc_term_rows),
            "run_csq_terms": pd.DataFrame(csq_term_rows),
        }

    rows = []
    mc_term_rows = []
    csq_term_rows = []

    for run_id, data in run_data.items():
        for row in data["variant_rows"]:
            mc_terms = _parse_mc_terms(row.get("mc", ""))
            csq_terms = _parse_csq_terms(row.get("csq", ""))
            for term in mc_terms:
                mc_term_rows.append({"run_id": run_id, "MC_term": term})
            for term in csq_terms:
                csq_term_rows.append({"run_id": run_id, "CSQ_term": term})

            rows.append(
                {
                    "run_id": run_id,
                    "gene_id": str(row.get("gene_id", "")),
                    "variant_key": row.get("variant_key", ""),
                    "mc_class": _classify_terms(mc_terms),
                    "csq_class": _classify_terms(csq_terms),
                }
            )

    df_run_consequence = pd.DataFrame(rows)
    if df_run_consequence.empty:
        empty = pd.DataFrame()
        return {
            "run_impact_class_shares": empty,
            "run_gene_impact_class_shares": empty,
            "run_mc_vs_csq_mismatch_summary": empty,
            "run_gene_mc_vs_csq_mismatch": empty,
            "run_mc_terms": empty,
            "run_csq_terms": empty,
        }

    df_mc_terms_by_run = (
        pd.DataFrame(mc_term_rows)
        .groupby(["run_id", "MC_term"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
        .sort_values(["run_id", "count"], ascending=[True, False])
        .reset_index(drop=True)
    ) if mc_term_rows else pd.DataFrame(columns=["run_id", "MC_term", "count"])

    df_csq_terms_by_run = (
        pd.DataFrame(csq_term_rows)
        .groupby(["run_id", "CSQ_term"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
        .sort_values(["run_id", "count"], ascending=[True, False])
        .reset_index(drop=True)
    ) if csq_term_rows else pd.DataFrame(columns=["run_id", "CSQ_term", "count"])

    long_rows = []
    for _, row in df_run_consequence.iterrows():
        if pd.notna(row["mc_class"]):
            long_rows.append(
                {"run_id": row["run_id"], "gene_id": row["gene_id"], "source": "ClinVar_MC", "impact_class": row["mc_class"]}
            )
        if pd.notna(row["csq_class"]):
            long_rows.append(
                {"run_id": row["run_id"], "gene_id": row["gene_id"], "source": "gnomAD_CSQ", "impact_class": row["csq_class"]}
            )

    df_impact_long = pd.DataFrame(long_rows)
    df_run_impact = (
        df_impact_long
        .groupby(["run_id", "source", "impact_class"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    df_run_source_totals = (
        df_run_impact
        .groupby(["run_id", "source"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "source_total"})
    )
    df_run_impact = df_run_impact.merge(df_run_source_totals, on=["run_id", "source"], how="left")
    df_run_impact["share_%"] = (df_run_impact["count"] / df_run_impact["source_total"] * 100).round(2)

    df_run_gene_impact = (
        df_impact_long
        .groupby(["run_id", "gene_id", "source", "impact_class"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    df_gene_source_totals = (
        df_run_gene_impact
        .groupby(["run_id", "gene_id", "source"], as_index=False)["count"]
        .sum()
        .rename(columns={"count": "source_total"})
    )
    df_run_gene_impact = df_run_gene_impact.merge(df_gene_source_totals, on=["run_id", "gene_id", "source"], how="left")
    df_run_gene_impact["share_%"] = (df_run_gene_impact["count"] / df_run_gene_impact["source_total"] * 100).round(2)

    df_both = df_run_consequence[df_run_consequence["mc_class"].notna() & df_run_consequence["csq_class"].notna()].copy()
    df_discordant = df_both[df_both["mc_class"] != df_both["csq_class"]].copy()

    if not df_both.empty:
        df_run_mismatch_summary = (
            df_discordant
            .groupby("run_id", as_index=False)
            .size()
            .rename(columns={"size": "discordant_count"})
        )
        df_both_summary = (
            df_both
            .groupby("run_id", as_index=False)
            .size()
            .rename(columns={"size": "both_annotated_count"})
        )
        df_run_mismatch_summary = df_run_mismatch_summary.merge(df_both_summary, on="run_id", how="right").fillna(0)
        df_run_mismatch_summary[["discordant_count", "both_annotated_count"]] = (
            df_run_mismatch_summary[["discordant_count", "both_annotated_count"]].astype(int)
        )
        df_run_mismatch_summary["discordant_share_%"] = (
            df_run_mismatch_summary["discordant_count"] / df_run_mismatch_summary["both_annotated_count"] * 100
        ).round(2)

        df_run_gene_mismatch = (
            df_discordant
            .groupby(["run_id", "gene_id"], as_index=False)
            .size()
            .rename(columns={"size": "discordant_count"})
        )
        df_run_gene_both = (
            df_both
            .groupby(["run_id", "gene_id"], as_index=False)
            .size()
            .rename(columns={"size": "both_annotated_count"})
        )
        df_run_gene_mismatch = df_run_gene_mismatch.merge(df_run_gene_both, on=["run_id", "gene_id"], how="left")
        df_run_gene_mismatch["discordant_share_%"] = (
            df_run_gene_mismatch["discordant_count"] / df_run_gene_mismatch["both_annotated_count"] * 100
        ).round(2)
        df_run_gene_mismatch = df_run_gene_mismatch.sort_values(["run_id", "discordant_count"], ascending=[True, False])
    else:
        df_run_mismatch_summary = pd.DataFrame()
        df_run_gene_mismatch = pd.DataFrame()

    return {
        "run_impact_class_shares": df_run_impact,
        "run_gene_impact_class_shares": df_run_gene_impact,
        "run_mc_vs_csq_mismatch_summary": df_run_mismatch_summary,
        "run_gene_mc_vs_csq_mismatch": df_run_gene_mismatch,
        "run_mc_terms": df_mc_terms_by_run,
        "run_csq_terms": df_csq_terms_by_run,
    }


def split_clnsig_labels(clnsig_raw: str) -> list[str]:
    if not clnsig_raw:
        return []
    return [item.strip() for item in re.split(r"[|,]", str(clnsig_raw)) if item.strip()]


def clinvar_stars_from_review_status(review_status: str) -> int | None:
    if not review_status:
        return None

    star_values = []
    for token in str(review_status).split("|"):
        token = token.strip().lower()
        if not token:
            continue
        if "practice_guideline" in token:
            star_values.append(4)
        elif "reviewed_by_expert_panel" in token:
            star_values.append(3)
        elif "criteria_provided,_multiple_submitters,_no_conflicts" in token:
            star_values.append(2)
        elif "criteria_provided,_single_submitter" in token or "criteria_provided,_conflicting_classifications" in token:
            star_values.append(1)
        elif (
            "no_assertion_criteria_provided" in token
            or "no_assertion_provided" in token
            or "no_interpretation_for_the_single_variant" in token
        ):
            star_values.append(0)
    return max(star_values) if star_values else None


def is_benign_label(label: str) -> bool:
    normalized = str(label).strip().lower()
    return ("benign" in normalized) and ("pathogenic" not in normalized) and ("conflicting" not in normalized)


def parse_clnsigconf_submitter_counts(clnsigconf_raw: str) -> defaultdict[str, int]:
    counts = defaultdict(int)
    if not clnsigconf_raw:
        return counts
    for token in str(clnsigconf_raw).split("|"):
        token = token.strip()
        if not token:
            continue
        match = re.match(r"^(.*)\((\d+)\)$", token)
        if match:
            label = match.group(1).strip()
            submitter_count = int(match.group(2))
        else:
            label = token
            submitter_count = 1
        if label:
            counts[label] += submitter_count
    return counts


def group_clnsig(clnsig: str) -> str:
    clnsig_lower = str(clnsig).lower()
    if "pathogenic" in clnsig_lower and "benign" not in clnsig_lower and "conflicting" not in clnsig_lower:
        return "Pathogenic"
    if "benign" in clnsig_lower and "pathogenic" not in clnsig_lower and "conflicting" not in clnsig_lower:
        return "Benign"
    if "uncertain" in clnsig_lower:
        return "VUS"
    if "conflicting" in clnsig_lower or "not_provided" in clnsig_lower:
        return "Conflicting"
    return "Other"


def _parse_mc_terms(raw: str) -> list[str]:
    out = []
    if not raw:
        return out
    for item in str(raw).split(","):
        item = item.strip()
        if not item:
            continue
        if "|" in item:
            term = item.split("|", 1)[1].strip()
        else:
            term = item
        if term:
            out.append(term)
    return out


def _parse_csq_terms(raw: str) -> list[str]:
    out = []
    if not raw:
        return out
    for item in re.split(r"[|,]", str(raw)):
        item = item.strip()
        if item:
            out.append(item)
    return out


def _classify_terms(terms: list[str]) -> str | None:
    if not terms:
        return None
    best = "OTHER"
    best_rank = IMPACT_RANK[best]
    for term in terms:
        impact = TERM_TO_IMPACT.get(term, "OTHER")
        rank = IMPACT_RANK[impact]
        if rank > best_rank:
            best = impact
            best_rank = rank
    return best


def _gene_sort_key(gene_id: str) -> tuple[int, str]:
    text = str(gene_id)
    return (0, f"{int(text):020d}") if text.isdigit() else (1, text)


def _summary_count_map(rows: list[dict], key_field: str) -> defaultdict[str, int]:
    counts = defaultdict(int)
    for row in rows:
        key = row.get(key_field)
        if key is None:
            continue
        counts[str(key)] += int(row.get("count", 0))
    return counts


if __name__ == "__main__":
    raise SystemExit(main())
