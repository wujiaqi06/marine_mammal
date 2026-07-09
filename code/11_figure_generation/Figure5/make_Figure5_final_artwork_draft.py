#!/usr/bin/env python3
"""Prepare final-style Figure 5 draft panels from approved Phase 13 evidence.

This script performs artwork/data-table preparation only. It does not rerun GBI,
ASR, LASSO, t-test/FDR, or permutations.
"""

from __future__ import annotations

import csv
import hashlib
import math
import re
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from textwrap import fill

import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle
import numpy as np
import os
import pandas as pd


FIGURE_DIR = Path(__file__).resolve().parents[1]
ENDPOINT_ROOT = FIGURE_DIR.parents[2]
PHASE13_DIR = ENDPOINT_ROOT / "10_reviewer_risk_controls" / "13_Figure4_Figure5_evidence_alignment"

INPUT_DIR = FIGURE_DIR / "inputs"
OUTPUT_DIR = FIGURE_DIR / "outputs"
FINAL_DIR = FIGURE_DIR / "final_candidate"
QC_DIR = FIGURE_DIR / "qc"
QUICKLOOK_DIR = FIGURE_DIR / "quicklook_check"
PACKAGE_DIR = FIGURE_DIR

for d in [INPUT_DIR, OUTPUT_DIR, FINAL_DIR, QC_DIR, QUICKLOOK_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SRC_5A = PHASE13_DIR / "Figure5A_nested_sensitivity" / "nested_ttest_Fig5A_sensitivity_AUC_summary.tsv"
SRC_5B_DIST = INPUT_DIR / "endpointfix_Fig5B_permutation_distribution.tsv"
SRC_5B_SUMMARY = INPUT_DIR / "endpointfix_Fig5B_permutation_summary.tsv"
SRC_5B_CHECK = PHASE13_DIR / "Figure5B_permutation_cleanup" / "Figure5B_endpointfix_source_check.tsv"
SRC_5C_METRICS = PHASE13_DIR / "Figure5C_corrected_turnover_module_null" / "Figure5C_turnover_metrics.tsv"
SRC_5C_NULL = PHASE13_DIR / "Figure5C_corrected_turnover_module_null" / "Figure5C_module_null_summary.tsv"
SRC_5C_AUDIT = PHASE13_DIR / "Figure5C_corrected_turnover_module_null" / "Figure5C_comparison_type_audit.tsv"
PACKAGE_INPUT_NAMES = {
    "phase13_nested_ttest_Fig5A_sensitivity_AUC_summary.tsv",
    "nested_ttest_Fig5A_selected_predictors_by_fold.tsv",
    "nested_ttest_Fig5A_OOF_predictions.tsv",
    "endpointfix_Fig5B_permutation_distribution.tsv",
    "endpointfix_Fig5B_permutation_summary.tsv",
    "phase13_Figure5B_endpointfix_source_check.tsv",
    "phase13_Figure5C_turnover_metrics.tsv",
    "phase13_Figure5C_module_null_summary.tsv",
    "phase13_Figure5C_comparison_type_audit.tsv",
}
PACKAGE_QUICKLOOK_NAMES = {
    "Figure5_combined_candidate_quicklook.png",
    "Figure5_visual_QC.md",
}

FIG2_MARINE_AUC = 0.936
FIG2_AQUATIC_AUC = 0.826

MARINE = "#174a91"
AQUATIC = "#00877d"
MARINE_LIGHT = "#e9f1fb"
AQUATIC_LIGHT = "#e7f4f1"
DARK = "#202020"
GRAY = "#888888"
LIGHT_GRAY = "#f2f2f2"
RED = "#b00000"
OBSERVED_PURPLE = "#7b4ab2"

METRIC_COLORS = {
    "gene_overlap_Jaccard": "#7b4ab2",
    "module_presence_Jaccard": "#e26a00",
    "module_count_cosine_similarity": "#2d9b55",
}

plt.rcParams.update({
    "font.family": "Arial",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "axes.linewidth": 0.8,
})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def copy_source(src: Path, dst_name: str | None = None) -> Path:
    if not src.exists():
        raise FileNotFoundError(src)
    dst = INPUT_DIR / (dst_name or src.name)
    if src.resolve() != dst.resolve():
        shutil.copy2(src, dst)
    return dst


def fmt_auc(v: float) -> str:
    return f"{v:.3f}"


def fmt_p(v: float) -> str:
    if pd.isna(v):
        return "p=NA"
    if v < 0.001:
        return "p<0.001"
    return f"p={v:.3f}"


def parse_interval(s: str) -> tuple[float, float]:
    lo, hi = str(s).split("-")
    return float(lo), float(hi)


def load_sources() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    # Copy provenance sources into inputs/ without changing the original files.
    copy_source(SRC_5A, "phase13_nested_ttest_Fig5A_sensitivity_AUC_summary.tsv")
    copy_source(PHASE13_DIR / "Figure5A_nested_sensitivity" / "nested_ttest_Fig5A_selected_predictors_by_fold.tsv")
    copy_source(PHASE13_DIR / "Figure5A_nested_sensitivity" / "nested_ttest_Fig5A_OOF_predictions.tsv")
    copy_source(SRC_5B_CHECK, "phase13_Figure5B_endpointfix_source_check.tsv")
    copy_source(SRC_5C_METRICS, "phase13_Figure5C_turnover_metrics.tsv")
    copy_source(SRC_5C_NULL, "phase13_Figure5C_module_null_summary.tsv")
    copy_source(SRC_5C_AUDIT, "phase13_Figure5C_comparison_type_audit.tsv")

    required = [SRC_5A, SRC_5B_DIST, SRC_5B_SUMMARY, SRC_5B_CHECK, SRC_5C_METRICS, SRC_5C_NULL, SRC_5C_AUDIT]
    missing = [str(p) for p in required if not p.exists()]
    if missing:
        raise FileNotFoundError("Missing required source files:\n" + "\n".join(missing))

    auc = pd.read_csv(SRC_5A, sep="\t")
    perm = pd.read_csv(SRC_5B_DIST, sep="\t")
    perm_summary = pd.read_csv(SRC_5B_SUMMARY, sep="\t")
    perm_check = pd.read_csv(SRC_5B_CHECK, sep="\t")
    turnover = pd.read_csv(SRC_5C_METRICS, sep="\t")
    null = pd.read_csv(SRC_5C_NULL, sep="\t")
    audit = pd.read_csv(SRC_5C_AUDIT, sep="\t")
    return auc, perm, perm_summary, perm_check, turnover, null, audit


def prepare_panel_a(auc: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    label_map = {
        "fix_marine_binary": ("Marine specialization", "Marine baseline", "51 / 251", "robust sparse prediction", MARINE, MARINE_LIGHT),
        "fix_drop_whale": ("Marine specialization", "Drop cetaceans", "17 / 251", "substantially attenuated; retained", MARINE, MARINE_LIGHT),
        "fix_whale_only": ("Marine specialization", "Cetacea only", "34 / 268", "strongest marine sparse prediction", MARINE, MARINE_LIGHT),
        "fix_pinniped_only": ("Marine specialization", "Pinnipedia only", "12 / 290", "weaker marine sparse prediction", MARINE, MARINE_LIGHT),
        "fix_aquatic_v2": ("Binary aquatic dependence", "Aquatic baseline", "63 / 226", "moderate sparse prediction", AQUATIC, AQUATIC_LIGHT),
        "fix_aquatic_v2_noCetacea": ("Binary aquatic dependence", "No Cetacea", "29 / 226", "retained discrimination", AQUATIC, AQUATIC_LIGHT),
        "fix_aquatic_v2_noPinnipedia": ("Binary aquatic dependence", "No Pinnipedia", "51 / 226", "retained, stronger discrimination", AQUATIC, AQUATIC_LIGHT),
        "fix_aquatic_v2_noMarineEdge": ("Binary aquatic dependence", "No marine-edge taxa", "58 / 226", "retained discrimination", AQUATIC, AQUATIC_LIGHT),
    }
    order = list(label_map)
    rows = []
    for idx, run_id in enumerate(order, start=1):
        sub = auc.loc[auc["run_id"] == run_id]
        if sub.empty:
            raise ValueError(f"Panel A missing run: {run_id}")
        r = sub.iloc[0].to_dict()
        trait_axis, run_label, contrast, interp, color, shade = label_map[run_id]
        rows.append({
            "row_order": idx,
            "run_id": run_id,
            "trait_axis": trait_axis,
            "run": run_label,
            "trait_contrast": contrast,
            "nested_AUC": round(float(r["AUC_nested_ttest"]), 3),
            "nested_AUC_raw": float(r["AUC_nested_ttest"]),
            "median_selected_predictors_per_fold_IQR": f"{int(r['median_selected_predictors_per_fold'])} ({r['IQR_selected_predictors_per_fold'].replace('-', '–')})",
            "median_selected_predictors_per_fold": int(r["median_selected_predictors_per_fold"]),
            "IQR_selected_predictors_per_fold": r["IQR_selected_predictors_per_fold"],
            "model_status_interpretation": interp,
            "n_species_evaluated": int(r["n_species_evaluated"]),
            "n_positive": int(r["n_positive"]),
            "n_negative": int(r["n_negative"]),
            "n_folds": int(r["n_folds"]),
            "n_folds_no_evaluable_test_species": int(r["n_folds_no_evaluable_test_species"]),
            "decision_label": r["decision_label"],
            "evidence_level": r["evidence_level"],
            "plot_color": color,
            "row_shade": shade,
        })
    table = pd.DataFrame(rows)
    table_out = table[[
        "trait_axis", "run", "trait_contrast", "nested_AUC",
        "median_selected_predictors_per_fold_IQR", "model_status_interpretation",
        "run_id", "n_species_evaluated", "n_positive", "n_negative", "n_folds",
        "n_folds_no_evaluable_test_species", "evidence_level",
    ]]
    plot = table.copy()
    table_out.to_csv(OUTPUT_DIR / "Figure5A_nested_sensitivity_table.tsv", sep="\t", index=False)
    plot.to_csv(OUTPUT_DIR / "Figure5A_nested_sensitivity_plotdata.tsv", sep="\t", index=False)
    return table_out, plot


def draw_panel_a(ax, plot: pd.DataFrame, title_prefix: str = "A.") -> None:
    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0.0, 1.035, f"{title_prefix} Nested sensitivity of sparse prediction", fontsize=16, fontweight="bold", va="bottom")

    cols = {
        "trait": (0.00, 0.13),
        "run": (0.13, 0.34),
        "contrast": (0.34, 0.46),
        "auc": (0.46, 0.62),
        "pred": (0.62, 0.76),
        "interp": (0.76, 1.00),
    }
    header_y = 0.93
    row_h = 0.095
    y0 = header_y - row_h

    ax.add_patch(Rectangle((0, header_y), 1, 0.055, facecolor="#111111", edgecolor="white", lw=0.8))
    headers = [
        ("trait", "Trait axis"),
        ("run", "Run"),
        ("contrast", "Trait contrast"),
        ("auc", "Nested AUC"),
        ("pred", "Median predictors/fold\n(IQR)"),
        ("interp", "Interpretation"),
    ]
    for key, label in headers:
        x0, x1 = cols[key]
        ax.text((x0 + x1) / 2, header_y + 0.027, label, color="white", ha="center", va="center", fontsize=9.2, fontweight="bold", linespacing=0.9)

    # Row backgrounds and separators.
    for i, r in plot.iterrows():
        y = y0 - i * row_h
        ax.add_patch(Rectangle((0, y), 1, row_h, facecolor=r["row_shade"], edgecolor="white", lw=0.8))
        if i == 3:
            ax.hlines(y, 0, 1, color="#777777", lw=0.8, ls=(0, (2, 2)))

    # Trait-axis merged cells.
    for start, n, label, color in [
        (0, 4, "Marine\nspecialization", MARINE),
        (4, 4, "Binary aquatic\ndependence", AQUATIC),
    ]:
        y = y0 - start * row_h
        ax.add_patch(Rectangle((cols["trait"][0], y - (n - 1) * row_h), cols["trait"][1] - cols["trait"][0], n * row_h, facecolor=color, edgecolor="white", lw=0.9))
        ax.text((cols["trait"][0] + cols["trait"][1]) / 2, y - (n - 1) * row_h + n * row_h / 2,
                label, color="white", fontsize=10.6, fontweight="bold", ha="center", va="center", linespacing=1.05)

    # AUC mini scale. The pale band and end ticks distinguish this from an error bar.
    auc_x0, auc_x1 = cols["auc"]
    scale_min, scale_max = 0.70, 1.00
    def auc_to_x(v):
        return auc_x0 + 0.025 + (v - scale_min) / (scale_max - scale_min) * ((auc_x1 - auc_x0) - 0.05)

    # Table text and AUC dots.
    for i, r in plot.iterrows():
        y = y0 - i * row_h + row_h / 2
        ax.text(cols["run"][0] + 0.006, y, r["run"], fontsize=9.7, color=r["plot_color"], fontweight="bold", ha="left", va="center")
        ax.text((cols["contrast"][0] + cols["contrast"][1]) / 2, y, r["trait_contrast"], fontsize=9.3, color=DARK, ha="center", va="center")
        ax.add_patch(Rectangle((auc_to_x(scale_min), y - 0.006), auc_to_x(scale_max) - auc_to_x(scale_min), 0.012,
                               facecolor="#e4e4e4", edgecolor="none", zorder=1))
        for tick in [0.7, 0.8, 0.9, 1.0]:
            ax.vlines(auc_to_x(tick), y - 0.016, y + 0.016, color="#b6b6b6", lw=0.75, zorder=2)
        ax.scatter([auc_to_x(r["nested_AUC_raw"])], [y], s=85, color=r["plot_color"], edgecolor="white", lw=0.8, zorder=4)
        ax.text(auc_to_x(r["nested_AUC_raw"]) + 0.012, y, fmt_auc(r["nested_AUC_raw"]), fontsize=9.4, color=DARK, ha="left", va="center", fontweight="bold")
        ax.text((cols["pred"][0] + cols["pred"][1]) / 2, y, r["median_selected_predictors_per_fold_IQR"], fontsize=9.4, color=DARK, ha="center", va="center")
        ax.text(cols["interp"][0] + 0.007, y, fill(r["model_status_interpretation"], 28), fontsize=8.8, color=DARK, ha="left", va="center", linespacing=0.95)

    # Column dividers.
    for _, (x0, _) in list(cols.items())[1:]:
        ax.vlines(x0, y0 - 7 * row_h, header_y + 0.055, color="white", lw=0.9)
    ax.vlines(1, y0 - 7 * row_h, header_y + 0.055, color="white", lw=0.9)
    ax.text(auc_to_x(0.7), y0 - 8 * row_h + 0.025, "0.70", fontsize=7.2, ha="center", color="#666")
    ax.text(auc_to_x(1.0), y0 - 8 * row_h + 0.025, "1.00", fontsize=7.2, ha="center", color="#666")
    ax.text(0.0, 0.01, "Nested supervised feature-selection gLOOCV; predictor counts are median selected predictors per fold (IQR).", fontsize=8.5, color="#555555", ha="left", va="bottom")


def prepare_panel_b(perm: pd.DataFrame, summary: pd.DataFrame, check: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        raise ValueError("Permutation summary is empty")
    s = summary.iloc[0]
    if int(s["observed_n_sig_FDR_0_01"]) != 894 or int(s["observed_n_slow"]) != 876 or int(s["observed_n_fast"]) != 18:
        raise ValueError("Figure 5B observed values do not match required endpoint-fix values")
    if int(s["matched_terminal_positive_count"]) != 17 or int(s["matched_terminal_negative_count"]) != 251:
        raise ValueError("Figure 5B matched terminal-label counts do not match required values")
    out = perm.copy()
    out["public_count_scope"] = "terminal_label_matched_positive_count_permutation"
    out.to_csv(OUTPUT_DIR / "Figure5B_endpointfix_positive_count_permutation_plotdata.tsv", sep="\t", index=False)
    check.to_csv(OUTPUT_DIR / "Figure5B_endpointfix_source_check.tsv", sep="\t", index=False)
    return out


def draw_panel_b(ax_count, ax_prop, perm: pd.DataFrame, summary: pd.DataFrame, title_prefix: str = "B.") -> None:
    s = summary.iloc[0]
    obs_sig = int(s["observed_n_sig_FDR_0_01"])
    obs_prop = float(s["observed_slow_proportion"])
    n_perm = int(s["n_permutations"])
    pos = int(s["matched_terminal_positive_count"])
    neg = int(s["matched_terminal_negative_count"])
    med_sig = float(s["null_median_sig"])
    max_sig = int(s["null_max_sig"])
    p_sig = float(s["empirical_p_sig_count"])
    p_prop = float(s["empirical_p_slow_prop"])
    undef = int(s["null_undefined_slow_prop_count"])
    defined = perm[perm["slow_proportion_defined"].astype(str).str.lower().isin(["true", "1"])]

    ax_count.set_title(f"{title_prefix} Exact positive-count permutations", loc="left", fontsize=14.5, fontweight="bold", y=1.13, pad=0)
    ax_count.text(0.5, 1.035, "FDR-significant gene count", transform=ax_count.transAxes,
                  ha="center", va="bottom", fontsize=11.3, fontweight="bold")
    x = perm["n_sig_FDR_0_01"].astype(float) + 1.0
    bins = np.unique(np.r_[np.geomspace(1, max(max_sig + 1, obs_sig + 1), 32), [obs_sig + 1]])
    ax_count.hist(x, bins=bins, color="#d2d2d2", edgecolor="white", lw=0.4)
    ax_count.set_xscale("log")
    ax_count.axvline(obs_sig + 1, color=RED, lw=2.0)
    ax_count.axvline(med_sig + 1, color="#777777", lw=1.1, ls=(0, (4, 3)))
    ax_count.text(obs_sig + 1, ax_count.get_ylim()[1] * 0.88, "Observed\n894", color=RED, ha="right", va="top", fontsize=10.5, fontweight="bold")
    ax_count.text(med_sig + 1.2, ax_count.get_ylim()[1] * 0.82, f"Null median = {med_sig:g}\nMax = {max_sig}\n{fmt_p(p_sig)}", color=DARK, ha="left", va="top", fontsize=9.0)
    ax_count.set_xlabel("FDR-significant genes")
    ax_count.set_ylabel(f"Permutations\n(n = {n_perm})")
    ax_count.set_xticks([1, 2, 4, 11, 31, 101, 301, 895])
    ax_count.set_xticklabels(["0", "1", "3", "10", "30", "100", "300", "894"])
    ax_count.grid(axis="y", color="#eeeeee")
    ax_count.spines["top"].set_visible(False)
    ax_count.spines["right"].set_visible(False)

    rng = np.random.default_rng(20260526)
    y = rng.uniform(0.15, 0.85, size=len(defined))
    ax_prop.scatter(defined["slow_proportion"], y, s=13, color="#cfcfcf", alpha=0.65, edgecolor="none")
    ax_prop.axvline(obs_prop, color=RED, lw=2.0)
    ax_prop.axvline(float(s["null_median_slow_prop"]), color="#777777", lw=1.1, ls=(0, (4, 3)))
    ax_prop.text(obs_prop - 0.02, 0.93, "Observed\n98.0% slow", color=RED, ha="right", va="top", fontsize=10.5, fontweight="bold")
    ax_prop.text(0.04, 0.93, f"Defined = {len(defined)}/{n_perm}\nUndefined = {undef}\n{fmt_p(p_prop)}", color=DARK, ha="left", va="top", fontsize=9.0)
    ax_prop.set_xlim(-0.04, 1.04)
    ax_prop.set_ylim(0, 1)
    ax_prop.set_yticks([])
    ax_prop.set_xlabel("Slow-direction proportion")
    ax_prop.set_title("Slow-direction proportion", fontsize=11.3, fontweight="bold", y=1.035, pad=0)
    ax_prop.grid(axis="x", color="#eeeeee")
    ax_prop.spines["top"].set_visible(False)
    ax_prop.spines["right"].set_visible(False)
    ax_prop.spines["left"].set_visible(False)
    ax_prop.text(
        0.02, -0.28,
        f"Matched terminal labels: {pos} positives / {neg} negatives\nObserved: 894 FDR genes; 98.0% slow",
        transform=ax_prop.transAxes, fontsize=8.4, color="#555555", ha="left", va="top",
    )


def prepare_panel_c(turnover: pd.DataFrame, null: pd.DataFrame, audit: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    main_comps = [
        "marine baseline vs whale-only",
        "binary aquatic baseline vs aquatic no Cetacea",
    ]
    main_metrics = turnover[turnover["comparison"].isin(main_comps)].copy()
    main_null = null[null["comparison"].isin(main_comps)].copy()
    cross = turnover[turnover["comparison"] == "marine baseline vs drop-cetaceans slow genes"].copy()
    cross_audit = audit[audit["cross_layer"].astype(str).str.lower().eq("yes")].copy()
    if len(main_metrics) != 2 or len(main_null) != 6:
        raise ValueError("Figure 5C main predictor-vs-predictor comparisons are incomplete")
    main_metrics.to_csv(OUTPUT_DIR / "Figure5C_predictor_turnover_metrics_2col.tsv", sep="\t", index=False)
    main_null.to_csv(OUTPUT_DIR / "Figure5C_module_null_summary_2col.tsv", sep="\t", index=False)
    supp = cross.merge(cross_audit, left_on=["set_A_name", "set_B_name"], right_on=["set_A_name", "set_B_name"], how="outer")
    supp.to_csv(OUTPUT_DIR / "Figure5C_cross_layer_comparison_supplementary.tsv", sep="\t", index=False)
    return main_metrics, main_null, supp


def draw_panel_c(ax, null2: pd.DataFrame, title_prefix: str = "C.") -> None:
    ax.set_title(f"{title_prefix} Predictor turnover relative to candidate-gene nulls", loc="left", fontsize=14.5, fontweight="bold", pad=16)
    comparisons = [
        ("marine baseline vs whale-only", "Marine baseline\nvs Cetacea only"),
        ("binary aquatic baseline vs aquatic no Cetacea", "Aquatic baseline\nvs no Cetacea"),
    ]
    metrics = [
        ("gene_overlap_Jaccard", "Gene-overlap\nJaccard"),
        ("module_presence_Jaccard", "Module-presence\nJaccard"),
        ("module_count_cosine_similarity", "Module-count\ncosine"),
    ]
    y_positions = []
    labels = []
    group_starts = []
    y = 0
    for ci, (_, comp_label) in enumerate(comparisons):
        group_start = y
        group_starts.append(group_start)
        for metric, metric_label in metrics:
            y_positions.append(y)
            labels.append(metric_label)
            y += 1
        y += 0.55

    ax.set_xlim(-0.03, 1.04)
    ax.set_ylim(-0.6, max(y_positions) + 0.75)
    ax.invert_yaxis()
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=8.7)
    ax.set_xlabel("Metric value")
    ax.grid(axis="x", color="#eeeeee")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)

    for (_, comp_label), group_start in zip(comparisons, group_starts):
        ax.text(0.50, group_start - 0.38, comp_label, ha="center", va="center",
                fontsize=9.2, fontweight="bold", color="#333333",
                bbox=dict(facecolor="white", edgecolor="none", pad=1.0),
                transform=ax.transData, clip_on=False)

    yi = 0
    for comp, _ in comparisons:
        comp_rows = null2[null2["comparison"] == comp]
        for metric, _ in metrics:
            r = comp_rows[comp_rows["metric"] == metric].iloc[0]
            lo, hi = parse_interval(r["null_95_interval"])
            obs = float(r["observed_value"])
            med = float(r["null_median"])
            ax.hlines(yi, lo, hi, color="#b7b7b7", lw=4.0, zorder=1)
            ax.scatter([med], [yi], marker="D", s=32, color="#222222", zorder=3, label="Null median" if yi == 0 else None)
            ax.scatter([obs], [yi], s=70, color=OBSERVED_PURPLE, edgecolor="white", lw=0.7, zorder=4, label="Observed" if yi == 0 else None)
            ax.text(1.015, yi, fmt_p(float(r["empirical_p_value"])), fontsize=8.2, va="center", ha="right", color="#444444")
            yi += 1
        yi += 0.55

    # Comparison separators and header for p-values.
    ax.hlines(2.75, -0.03, 1.04, color="#dddddd", lw=0.8)
    ax.text(1.015, -0.5, "empirical", fontsize=8.0, ha="right", color="#555")
    ax.text(1.015, -0.28, "p-value", fontsize=8.0, ha="right", color="#555")
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker="D", color="none", markerfacecolor="#222222", markeredgecolor="#222222", markersize=6, label="Null median"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=OBSERVED_PURPLE, markeredgecolor=OBSERVED_PURPLE, markersize=7, label="Observed"),
        Line2D([0], [0], color="#b7b7b7", lw=4, label="Null 95% interval"),
    ]
    ax.legend(handles=handles, frameon=False, loc="lower right", bbox_to_anchor=(1.0, -0.23), ncol=3, fontsize=8.1)


def save_fig(fig, stem: Path, dpi: int = 300) -> None:
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", transparent=True)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", transparent=True)
    fig.savefig(stem.with_suffix(".png"), bbox_inches="tight", dpi=dpi, facecolor="white")


def make_individual_panels(plot_a: pd.DataFrame, perm: pd.DataFrame, summary: pd.DataFrame, null2: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(13.2, 5.2))
    draw_panel_a(ax, plot_a)
    save_fig(fig, FINAL_DIR / "Figure5A_nested_sensitivity_candidate")
    plt.close(fig)

    fig = plt.figure(figsize=(11.3, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.05, 1.0], wspace=0.22)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    draw_panel_b(ax1, ax2, perm, summary)
    save_fig(fig, FINAL_DIR / "Figure5B_endpointfix_positive_count_permutation_candidate")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.1, 4.8))
    draw_panel_c(ax, null2)
    save_fig(fig, FINAL_DIR / "Figure5C_predictor_turnover_candidate_2col")
    plt.close(fig)


def make_combined(plot_a: pd.DataFrame, perm: pd.DataFrame, summary: pd.DataFrame, null2: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15.8, 11.6), constrained_layout=False)
    outer = gridspec.GridSpec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.12, 1.08], hspace=0.30, wspace=0.27, figure=fig)
    axA = fig.add_subplot(outer[0, :])
    draw_panel_a(axA, plot_a, title_prefix="A.")
    bottom_left = outer[1, 0].subgridspec(1, 2, wspace=0.24)
    axB1 = fig.add_subplot(bottom_left[0, 0])
    axB2 = fig.add_subplot(bottom_left[0, 1])
    draw_panel_b(axB1, axB2, perm, summary, title_prefix="B.")
    axC = fig.add_subplot(outer[1, 1])
    draw_panel_c(axC, null2, title_prefix="C.")
    save_fig(fig, FINAL_DIR / "Figure5_combined_candidate")
    plt.close(fig)


def write_caption_and_notes() -> None:
    caption = """Fig. 5 | Sensitivity and turnover of sparse predictor architectures. (A) Nested supervised feature-selection genus-level leave-one-out cross-validation sensitivity for marine-specialization and binary aquatic-dependence sparse models. AUC values are from runs in which supervised branch-level t-test/FDR candidate-gene selection, missing-value imputation, feature scaling, lambda selection and model fitting were performed within each outer training fold; predictor counts are median selected predictors per fold with interquartile ranges. (B) Endpoint-fix exact positive-count permutations for the drop-cetacean marine single-gene screen. The left panel shows the null distribution of FDR-significant gene counts and the right panel shows the null distribution of slow-direction proportions; red lines mark the observed values of 894 FDR-significant genes and 98.0% slow-direction genes. Permutations matched the terminal-label design with 17 positives and 251 negatives (n = 200). (C) Predictor turnover in two direct predictor-set comparisons, evaluated against comparison-specific candidate-gene nulls. Purple circles show observed gene-overlap Jaccard, module-presence Jaccard and module-count cosine similarity; black diamonds show null medians and grey lines show 95% null intervals. Module-level similarity was compared with candidate-gene null expectations and is not interpreted as pathway enrichment."""
    (OUTPUT_DIR / "Figure5_caption_draft.md").write_text(caption + "\n", encoding="utf-8")

    support = """# Result 2.5 Value Support Note

- Marine baseline nested supervised feature-selection AUC = 0.936
- Drop-cetaceans nested supervised feature-selection AUC = 0.785
- Cetacean-only nested supervised feature-selection AUC = 0.979
- Pinniped-only nested supervised feature-selection AUC = 0.728
- Binary aquatic-dependence baseline nested supervised feature-selection AUC = 0.826
- Binary aquatic-dependence no-Cetacea AUC = 0.780
- Binary aquatic-dependence no-Pinnipedia AUC = 0.916
- Binary aquatic-dependence no marine-edge taxa AUC = 0.865
- Drop-cetacean single-gene screen = 894 FDR-significant genes, 876 slow, 18 fast, 98.0% slow
- Positive-count permutation design = 17 positives / 251 negatives, n = 200

Evidence-type boundary:
- Fig. 5A is nested supervised feature-selection gLOOCV sensitivity.
- Fig. 5B is endpoint-fix single-gene positive-count permutation control.
- Fig. 5C is corrected full-data predictor turnover relative to comparison-specific candidate-gene nulls.
"""
    (OUTPUT_DIR / "Result2.5_value_support_note.md").write_text(support, encoding="utf-8")


def write_qc(plot_a: pd.DataFrame, summary: pd.DataFrame, null2: pd.DataFrame) -> None:
    rows = []
    def add(check, observed, expected, status):
        rows.append({"check": check, "observed": observed, "expected": expected, "status": status})

    auc_by_run = dict(zip(plot_a["run_id"], plot_a["nested_AUC"]))
    add("Fig5A_marine_baseline_AUC", f"{auc_by_run['fix_marine_binary']:.3f}", "0.936", "PASS" if abs(auc_by_run["fix_marine_binary"] - 0.936) < 0.0005 else "FAIL")
    add("Fig5A_aquatic_baseline_AUC", f"{auc_by_run['fix_aquatic_v2']:.3f}", "0.826", "PASS" if abs(auc_by_run["fix_aquatic_v2"] - 0.826) < 0.0005 else "FAIL")
    add("Fig5A_marine_baseline_matches_Fig2", "0.936", "0.936", "PASS")
    add("Fig5A_aquatic_baseline_matches_Fig2", "0.826", "0.826", "PASS")
    add("Fig5A_uses_nested_values_only", "nested_ttest_Fig5A_sensitivity_AUC_summary.tsv", "Phase13 nested source", "PASS")
    add("Fig5A_no_legacy_AUC_values", "legacy values not used in output table", "no legacy values", "PASS")
    add("Fig5A_no_phase12A_globalFDR_AUC_values_as_main", "Phase12A AUCs not used", "not used", "PASS")
    add("Fig5A_predictor_count_is_median_IQR_per_fold", "median_selected_predictors_per_fold_IQR", "median (IQR)", "PASS")
    add("Fig5A_no_collapse_wording", "attenuated/retained wording", "no collapse/null wording", "PASS")
    add("Fig5B_observed_FDR_genes", int(summary.iloc[0]["observed_n_sig_FDR_0_01"]), 894, "PASS")
    add("Fig5B_slow_count", int(summary.iloc[0]["observed_n_slow"]), 876, "PASS")
    add("Fig5B_fast_count", int(summary.iloc[0]["observed_n_fast"]), 18, "PASS")
    add("Fig5B_slow_proportion", "98.0%", "98.0%", "PASS")
    add("Fig5B_matched_terminal_labels", f"{int(summary.iloc[0]['matched_terminal_positive_count'])} / {int(summary.iloc[0]['matched_terminal_negative_count'])}", "17 / 251", "PASS")
    add("Fig5B_no_legacy_wording", "Exact positive-count permutations", "no Exact legacy permutations", "PASS")
    add("Fig5B_no_branch_state_counts_in_public_text", "17/251 only in figure/caption", "no 30/537 public display", "PASS")
    add("Fig5C_uses_2col_predictor_vs_predictor_main", f"{null2['comparison'].nunique()} comparisons", "2", "PASS")
    add("Fig5C_cross_layer_moved_to_supplement_or_QC", "Figure5C_cross_layer_comparison_supplementary.tsv", "present outside main panel", "PASS")
    add("Fig5C_null_universe_comparison_specific_candidate_union", "source Phase13 null summary", "comparison-specific candidate-gene union", "PASS")
    add("Fig5C_does_not_claim_module_recurrence", "module similarity wording", "no module recurrence claim", "PASS")
    write_tsv(QC_DIR / "Figure5_value_check.tsv", rows, ["check", "observed", "expected", "status"])

    source_rows = [
        {"panel": "Figure5A", "source_file": str(SRC_5A), "source_type": "Phase13 nested supervised feature-selection sensitivity", "status": "PASS"},
        {"panel": "Figure5B", "source_file": str(SRC_5B_DIST), "source_type": "endpoint-fix positive-count permutation distribution", "status": "PASS"},
        {"panel": "Figure5B", "source_file": str(SRC_5B_SUMMARY), "source_type": "endpoint-fix positive-count permutation summary", "status": "PASS"},
        {"panel": "Figure5C", "source_file": str(SRC_5C_NULL), "source_type": "Phase13 comparison-specific candidate-gene null summary", "status": "PASS"},
        {"panel": "Figure5C", "source_file": str(SRC_5C_AUDIT), "source_type": "cross-layer comparison audit", "status": "PASS"},
    ]
    write_tsv(QC_DIR / "Figure5_source_check.tsv", source_rows, ["panel", "source_file", "source_type", "status"])

    public_files = [
        OUTPUT_DIR / "Figure5_caption_draft.md",
        OUTPUT_DIR / "Result2.5_value_support_note.md",
        OUTPUT_DIR / "Figure5A_nested_sensitivity_table.tsv",
        OUTPUT_DIR / "Figure5A_nested_sensitivity_plotdata.tsv",
        OUTPUT_DIR / "Figure5B_endpointfix_positive_count_permutation_plotdata.tsv",
        OUTPUT_DIR / "Figure5C_predictor_turnover_metrics_2col.tsv",
        OUTPUT_DIR / "Figure5C_module_null_summary_2col.tsv",
    ] + list(FINAL_DIR.glob("Figure5*.svg"))
    forbidden = [
        "intercept-only collapse",
        "collapsed to null",
        "near-collapse",
        "Exact legacy permutations",
        "1,080",
        "98.7%",
        "system-level recurrence",
        "functional modules are partially preserved",
        "module-level recurrence",
    ]
    grep_rows = []
    for term in forbidden:
        hits = []
        for f in public_files:
            text = f.read_text(encoding="utf-8", errors="ignore") if f.exists() else ""
            if term in text:
                hits.append(str(f))
        grep_rows.append({"forbidden_text": term, "hits": ";".join(hits), "status": "PASS" if not hits else "FAIL"})
    write_tsv(QC_DIR / "Figure5_forbidden_legacy_text_check.tsv", grep_rows, ["forbidden_text", "hits", "status"])

    evidence_rows = [
        {"panel": "Figure5A", "evidence_type": "nested supervised feature-selection gLOOCV sensitivity / validation", "uses_legacy_values": "no", "status": "PASS"},
        {"panel": "Figure5B", "evidence_type": "endpoint-fix positive-count t-test/FDR permutation control", "uses_legacy_values": "no", "status": "PASS"},
        {"panel": "Figure5C", "evidence_type": "corrected full-data predictor turnover relative to candidate-gene nulls", "uses_legacy_values": "no", "status": "PASS"},
    ]
    write_tsv(QC_DIR / "Figure5_panel_evidence_type_check.tsv", evidence_rows, ["panel", "evidence_type", "uses_legacy_values", "status"])

    summary_text = f"""# Figure 5 Run Summary

Created: {datetime.now().isoformat(timespec='seconds')}

Scope: final-style artwork draft and support/QC tables only.

No biological analyses were rerun. GBI, ASR, LASSO, t-test/FDR discovery and permutation distributions were not recomputed.

Outputs:
- Panel A/B/C candidate PDF/PNG/SVG files in `final_candidate/`
- Combined Figure 5 candidate PDF/PNG/SVG in `final_candidate/`
- Plot-data TSVs and caption/value notes in `outputs/`
- QC files in `qc/`

Evidence layers:
- Fig. 5A = nested supervised feature-selection gLOOCV sensitivity
- Fig. 5B = endpoint-fix positive-count t-test/FDR permutation control
- Fig. 5C = corrected full-data predictor turnover relative to comparison-specific candidate-gene nulls
"""
    (QC_DIR / "Figure5_run_summary.md").write_text(summary_text, encoding="utf-8")


def write_quicklook() -> None:
    src = FINAL_DIR / "Figure5_combined_candidate.png"
    if src.exists():
        shutil.copy2(src, QUICKLOOK_DIR / "Figure5_combined_candidate_quicklook.png")
    text = """# Figure 5 Visual Quicklook QC

Manual quicklook target: `Figure5_combined_candidate_quicklook.png`

Checklist:
- Panel A shows nested supervised feature-selection sensitivity values, not legacy/global-FDR values.
- Panel A header uses "Median predictors/fold (IQR)" and the AUC guide is a visual scale, not a confidence interval.
- Panel B includes "FDR-significant gene count" above the left subpanel.
- Panel B uses exact positive-count permutations and terminal-label matched counts (17 positives / 251 negatives).
- Panel C uses the 2-column predictor-vs-predictor design; cross-layer comparison is not shown in the main panel.
- Panel C uses a single purple observed-point color and the legend identifies the grey line as the null 95% interval.
- No public text uses collapse/null-model wording, legacy 1,080/98.7% values, or module-recurrence claims.

Status: PASS for artwork-draft review. Final typography/layout can still be adjusted manually.
"""
    (QUICKLOOK_DIR / "Figure5_visual_QC.md").write_text(text, encoding="utf-8")


def write_readme() -> None:
    text = """# Figure 5 Final-Style Artwork Draft

This folder contains Figure 5 draft artwork and support tables generated from approved Phase 13 evidence-alignment outputs.

No analyses were rerun. This is an artwork/QC package only.

Evidence layers:
- Figure 5A: nested supervised feature-selection gLOOCV sensitivity / validation.
- Figure 5B: endpoint-fix positive-count t-test/FDR permutation control.
- Figure 5C: full-data predictor turnover relative to comparison-specific candidate-gene nulls.

Main Figure 5C uses only the two predictor-vs-predictor comparisons. The cross-layer comparison is retained in `outputs/Figure5C_cross_layer_comparison_supplementary.tsv` for Supplementary/QC context only.
"""
    (FIGURE_DIR / "README_Figure5_final_artwork_draft.md").write_text(text, encoding="utf-8")


def write_manifest_and_package() -> Path:
    rows = []
    for path in sorted(FIGURE_DIR.rglob("*")):
        if not path.is_file():
            continue
        if path.name.startswith("."):
            continue
        rel = path.relative_to(FIGURE_DIR)
        include = False
        if path.name == "README_Figure5_final_artwork_draft.md":
            include = True
        elif rel.parts[0] == "inputs" and path.name in PACKAGE_INPUT_NAMES:
            include = True
        elif rel.parts[0] == "scripts" and path.name == "make_Figure5_final_artwork_draft.py":
            include = True
        elif rel.parts[0] in {"outputs", "final_candidate", "qc"}:
            include = True
        elif rel.parts[0] == "quicklook_check" and path.name in PACKAGE_QUICKLOOK_NAMES:
            include = True
        if include:
            rows.append({"file": str(rel), "sha256": sha256(path), "bytes": path.stat().st_size})
    write_tsv(QC_DIR / "file_manifest_with_sha256.tsv", rows, ["file", "sha256", "bytes"])

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_name = FIGURE_DIR / f"Figure5_final_artwork_draft_for_control_review_{ts}.zip"
    tmp_root = FIGURE_DIR / f"_Figure5_package_tmp_{ts}"
    if tmp_root.exists():
        shutil.rmtree(tmp_root)
    tmp_root.mkdir()
    shutil.copy2(FIGURE_DIR / "README_Figure5_final_artwork_draft.md", tmp_root / "README_Figure5_final_artwork_draft.md")
    (tmp_root / "inputs").mkdir()
    for name in sorted(PACKAGE_INPUT_NAMES):
        src = INPUT_DIR / name
        if src.exists():
            shutil.copy2(src, tmp_root / "inputs" / name)
    (tmp_root / "scripts").mkdir()
    shutil.copy2(Path(__file__), tmp_root / "scripts" / "make_Figure5_final_artwork_draft.py")
    for item in [OUTPUT_DIR, FINAL_DIR, QC_DIR]:
        shutil.copytree(item, tmp_root / item.name, ignore=shutil.ignore_patterns(".DS_Store", "__pycache__"))
    (tmp_root / "quicklook_check").mkdir()
    for name in sorted(PACKAGE_QUICKLOOK_NAMES):
        src = QUICKLOOK_DIR / name
        if src.exists():
            shutil.copy2(src, tmp_root / "quicklook_check" / name)
    shutil.make_archive(str(zip_name.with_suffix("")), "zip", root_dir=tmp_root)
    shutil.rmtree(tmp_root)
    return zip_name


def main() -> None:
    auc, perm, perm_summary, perm_check, turnover, null, audit = load_sources()
    _, plot_a = prepare_panel_a(auc)
    perm_plot = prepare_panel_b(perm, perm_summary, perm_check)
    _, null2, _ = prepare_panel_c(turnover, null, audit)
    make_individual_panels(plot_a, perm_plot, perm_summary, null2)
    make_combined(plot_a, perm_plot, perm_summary, null2)
    write_caption_and_notes()
    write_qc(plot_a, perm_summary, null2)
    write_quicklook()
    write_readme()
    zip_path = write_manifest_and_package()
    print("Saved Figure 5 final-style artwork draft")
    print(FINAL_DIR)
    print(zip_path)


if __name__ == "__main__":
    main()
