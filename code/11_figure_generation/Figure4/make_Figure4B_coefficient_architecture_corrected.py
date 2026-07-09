#!/usr/bin/env python3
"""
Figure 4B: Coefficient architecture of LASSO predictors.

This is a custom matplotlib diagram rather than a standard statistical plot.
It uses baseline marine_binary and aquatic_v2 LASSO coefficients only.

Coefficient sign convention used for color:
- The LASSO model is fit to predict trait = 1 from GBI predictors.
- Positive coefficient: higher GBI increases predicted trait log-odds.
- Under the current Figure 4 visual convention, positive coefficients are
  shown as fast-direction / higher-rate predictors in red.
- Negative coefficients are shown as slow-direction / lower-rate predictors
  in blue.

No t-test values are used in this panel.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
from matplotlib.patches import Patch


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
INPUT_DIR = FIGURE_DIR / "inputs"
OUTPUT_DIR = FIGURE_DIR / "Figure4B_outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

GENE_TABLE_FILE = INPUT_DIR / "Table1_gene_level_working_table.tsv"
PLOT_TABLE_FILE = OUTPUT_DIR / "Figure4B_plot_table.tsv"
PDF_FILE = OUTPUT_DIR / "Figure4B_coefficient_architecture.pdf"
SVG_FILE = OUTPUT_DIR / "Figure4B_coefficient_architecture.svg"
PNG_FILE = OUTPUT_DIR / "Figure4B_coefficient_architecture.png"

BLUE = "#1F4E79"
RED = "#8B0000"
MARINE_BLUE = "#4C78A8"
AQUATIC_TEAL = "#44A08D"
SHARED_GREY = "#9E9E9E"
LIGHT_GREY = "#D6D6D6"
DARK = "#111111"

REPRESENTATIVE_GENES = [
    "TGM3", "KRT15", "ALOX12B", "A2ML1",
    "EPOR", "CD34", "F11", "ITGA2B",
    "TREM2", "HHLA1", "IL33", "MEFV",
    "OR52W1", "OMP", "PDE6C", "OPN1SW",
    "GPR68", "KCNK18", "CALHM1", "PIEZO1",
    "CATSPER4", "CATSPERZ", "CATSPERG",
    "TPX2", "FANCM", "KIF18A",
    "MYH13", "COL4A1", "COL20A1",
]


def choose_representative_genes(df: pd.DataFrame, n_each: int = 12) -> list[str]:
    """Select corrected-data labels rather than carrying legacy label choices."""
    labels: set[str] = set()
    for group, coef_col in [
        ("marine_specific", "abs_marine_coef"),
        ("aquatic_specific", "abs_aquatic_coef"),
        ("shared", "max_abs_coef"),
    ]:
        sub = df[df["lasso_group"].eq(group)].copy()
        if sub.empty:
            continue
        sub = sub.sort_values([coef_col, "gene"], ascending=[False, True])
        labels.update(sub.head(n_each)["gene"].tolist())
    return sorted(labels)


def require_columns(df: pd.DataFrame, cols: list[str], file_label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {file_label}: {', '.join(missing)}")


def as_bool(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().eq("TRUE")


def coef_color(coef: float) -> str:
    return RED if coef >= 0 else BLUE


def coef_direction(coef: float) -> str:
    return "fast_direction_higher_rate" if coef >= 0 else "slow_direction_lower_rate"


def point_size(abs_coef: float) -> float:
    if abs_coef > 0.5:
        return 42
    if abs_coef > 0.2:
        return 30
    if abs_coef > 0.1:
        return 22
    return 11


def spread_label_positions(y_values: pd.Series, min_gap: float = 2.2,
                           lower: float = 2.0, upper: float = 82.0) -> np.ndarray:
    """Simple deterministic vertical label spreading for side labels."""
    y = np.asarray(y_values, dtype=float)
    order = np.argsort(y)
    adjusted = y.copy()
    last = lower - min_gap
    for idx in order:
        adjusted[idx] = max(y[idx], last + min_gap)
        last = adjusted[idx]
    overflow = adjusted[order[-1]] - upper if len(order) else 0
    if overflow > 0:
        adjusted -= overflow
        # A second forward pass keeps the lower bound and spacing sensible.
        last = lower - min_gap
        for idx in order:
            adjusted[idx] = max(adjusted[idx], last + min_gap)
            last = adjusted[idx]
    return adjusted


def model_point_table(df: pd.DataFrame, model: str, axis_x: float, max_abs: float) -> pd.DataFrame:
    if model == "marine":
        selected_col = "selected_in_marine"
        coef_col = "marine_coef"
        abs_col = "abs_marine_coef"
        model_label = "Marine model"
    else:
        selected_col = "selected_in_aquatic"
        coef_col = "aquatic_coef"
        abs_col = "abs_aquatic_coef"
        model_label = "Aquatic model"

    sub = df[df[selected_col]].copy()
    sub[coef_col] = pd.to_numeric(sub[coef_col], errors="coerce")
    sub[abs_col] = pd.to_numeric(sub[abs_col], errors="coerce")
    sub = sub.dropna(subset=[coef_col])
    sub = sub[sub[coef_col] != 0].copy()

    # Positive coefficients are placed higher than negative coefficients, with
    # larger absolute coefficients farther from the middle.
    sub = sub.sort_values([coef_col, abs_col, "gene"], ascending=[False, False, True]).reset_index(drop=True)
    n = len(sub)
    sub["y"] = np.linspace(n, 1, n)
    sub["model"] = model_label
    sub["coef"] = sub[coef_col]
    sub["abs_coef"] = sub[abs_col]
    sub["axis_x"] = axis_x
    sub["x"] = axis_x + (sub["coef"] / max_abs) * 0.86
    sub["coefficient_direction"] = sub["coef"].map(coef_direction)
    sub["coefficient_sign"] = np.where(sub["coef"] >= 0, "positive", "negative")
    sub["color"] = sub["coef"].map(coef_color)
    sub["size"] = sub["abs_coef"].map(point_size)
    sub["is_representative_label"] = sub["gene"].isin(REPRESENTATIVE_GENES)
    sub["gene_set"] = np.where(
        sub["lasso_group"].eq("shared"),
        "Shared aquatic foundation",
        np.where(sub["lasso_group"].eq("marine_specific"), "Marine-specific predictors", "Aquatic-specific predictors"),
    )
    return sub[[
        "gene", "model", "lasso_group", "gene_set", "coef", "abs_coef", "coefficient_sign",
        "coefficient_direction", "axis_x", "x", "y", "size", "is_representative_label"
    ]]


def main() -> None:
    required = [
        "gene", "lasso_group", "marine_coef", "aquatic_coef",
        "abs_marine_coef", "abs_aquatic_coef", "selected_in_marine", "selected_in_aquatic",
    ]
    df = pd.read_csv(GENE_TABLE_FILE, sep="\t")
    require_columns(df, required, str(GENE_TABLE_FILE))
    df["selected_in_marine"] = as_bool(df["selected_in_marine"])
    df["selected_in_aquatic"] = as_bool(df["selected_in_aquatic"])
    for col in ["marine_coef", "aquatic_coef", "abs_marine_coef", "abs_aquatic_coef"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["max_abs_coef"] = df[["abs_marine_coef", "abs_aquatic_coef"]].max(axis=1, skipna=True)

    marine_n = int(df["selected_in_marine"].sum())
    aquatic_n = int(df["selected_in_aquatic"].sum())
    shared_n = int((df["lasso_group"] == "shared").sum())
    marine_only_n = int((df["lasso_group"] == "marine_specific").sum())
    aquatic_only_n = int((df["lasso_group"] == "aquatic_specific").sum())
    print("Figure 4B input summary")
    print(f"marine genes = {marine_n}")
    print(f"aquatic genes = {aquatic_n}")
    print(f"shared genes = {shared_n}")
    print(f"marine-specific genes = {marine_only_n}")
    print(f"aquatic-specific genes = {aquatic_only_n}")
    print("Coefficient sign convention: positive = red/higher-rate; negative = blue/lower-rate.")

    global REPRESENTATIVE_GENES
    REPRESENTATIVE_GENES = choose_representative_genes(df, n_each=12)

    max_abs = float(np.nanmax(np.abs(pd.concat([df["marine_coef"], df["aquatic_coef"]]))))
    if not math.isfinite(max_abs) or max_abs <= 0:
        raise ValueError("No nonzero LASSO coefficients found.")

    marine_axis_x = -2.85
    aquatic_axis_x = 2.85
    shared_x = 0.0

    marine_pts = model_point_table(df, "marine", marine_axis_x, max_abs)
    aquatic_pts = model_point_table(df, "aquatic", aquatic_axis_x, max_abs)

    shared = df[df["lasso_group"] == "shared"].copy()
    shared["shared_mean_coef"] = shared[["marine_coef", "aquatic_coef"]].mean(axis=1)
    shared["shared_abs_mean_coef"] = shared["shared_mean_coef"].abs()
    shared["shared_color"] = shared["shared_mean_coef"].map(coef_color)
    shared["shared_coefficient_direction"] = shared["shared_mean_coef"].map(coef_direction)
    shared["marine_aquatic_sign_opposite"] = (
        np.sign(shared["marine_coef"]) * np.sign(shared["aquatic_coef"]) < 0
    )
    n_opposite = int(shared["marine_aquatic_sign_opposite"].sum())
    if n_opposite:
        print(f"WARNING: {n_opposite} shared genes have opposite marine/aquatic coefficient signs.")
    shared = shared.sort_values(
        ["shared_mean_coef", "shared_abs_mean_coef", "gene"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    shared["shared_x"] = shared_x
    shared["shared_y"] = np.linspace(max(marine_n, aquatic_n) - 6, 8, len(shared))

    shared_lookup = shared[[
        "gene", "shared_x", "shared_y", "shared_mean_coef", "shared_abs_mean_coef",
        "shared_color", "shared_coefficient_direction", "marine_aquatic_sign_opposite"
    ]]
    marine_pts = marine_pts.merge(shared_lookup, on="gene", how="left")
    aquatic_pts = aquatic_pts.merge(shared_lookup, on="gene", how="left")
    plot_table = pd.concat([marine_pts, aquatic_pts], ignore_index=True)
    plot_table.to_csv(PLOT_TABLE_FILE, sep="\t", index=False)

    fig, ax = plt.subplots(figsize=(9.6, 6.2))
    ax.set_facecolor("white")

    # faint horizontal guides
    for y in np.linspace(5, max(marine_n, aquatic_n) - 2, 22):
        ax.hlines(y, marine_axis_x - 0.95, aquatic_axis_x + 0.95, color="#EFEFEF", lw=0.45, zorder=0)

    # model-specific zero axes
    ax.axvline(marine_axis_x, color="#444444", lw=0.8, ls=(0, (4, 4)), zorder=1)
    ax.axvline(aquatic_axis_x, color="#444444", lw=0.8, ls=(0, (4, 4)), zorder=1)

    # shared-gene center column and connectors
    for _, row in shared.iterrows():
        gene = row["gene"]
        cy = row["shared_y"]
        m = marine_pts[marine_pts["gene"].eq(gene)]
        a = aquatic_pts[aquatic_pts["gene"].eq(gene)]
        if len(m) == 1:
            ax.plot([m.iloc[0]["x"], shared_x - 0.12], [m.iloc[0]["y"], cy],
                    color=LIGHT_GREY, lw=0.7, zorder=1)
        if len(a) == 1:
            ax.plot([shared_x + 0.12, a.iloc[0]["x"]], [cy, a.iloc[0]["y"]],
                    color=LIGHT_GREY, lw=0.7, zorder=1)
        ax.scatter(shared_x, cy, s=18, color=SHARED_GREY, edgecolor="white", linewidth=0.3, zorder=3)

    # all model points
    for pts, marker in [(marine_pts, "^"), (aquatic_pts, "s")]:
        ax.scatter(
            pts["x"], pts["y"],
            s=pts["size"],
            c=pts["coef"].map(coef_color),
            marker=marker,
            alpha=np.where(pts["is_representative_label"], 0.98, 0.42),
            edgecolor="none",
            zorder=4,
        )

    # central labels for shared genes: label all shared names, small and aligned.
    for _, row in shared.iterrows():
        ax.text(shared_x + 0.07, row["shared_y"], row["gene"], va="center", ha="left",
                fontsize=7.0, color=row["shared_color"])

    # representative labels on side plots
    def add_side_labels(pts: pd.DataFrame, side: str) -> None:
        label_pts = pts[pts["is_representative_label"]].copy()
        if label_pts.empty:
            return
        label_pts = label_pts.sort_values("y").copy()
        label_pts["label_y"] = spread_label_positions(label_pts["y"], min_gap=2.35, lower=3.0, upper=max(marine_n, aquatic_n) + 1)
        for _, row in label_pts.iterrows():
            coef = row["coef"]
            if side == "marine":
                ha = "right" if coef < 0 else "left"
                dx = -0.06 if coef < 0 else 0.06
            else:
                ha = "left" if coef > 0 else "right"
                dx = 0.06 if coef > 0 else -0.06
            if abs(row["label_y"] - row["y"]) > 0.3:
                ax.plot([row["x"], row["x"] + dx * 0.8], [row["y"], row["label_y"]],
                        color="#BDBDBD", lw=0.45, zorder=2)
            ax.text(row["x"] + dx, row["label_y"], row["gene"], va="center", ha=ha,
                    fontsize=7.4, color=coef_color(coef),
                    fontstyle="italic" if len(row["gene"]) > 5 else "normal")

    add_side_labels(marine_pts, "marine")
    add_side_labels(aquatic_pts, "aquatic")

    # Titles and annotations
    top_y = max(marine_n, aquatic_n) + 5
    ax.text(marine_axis_x, top_y, f"Marine model\n(n = {marine_n})", ha="center", va="bottom",
            color=BLUE, fontsize=10.5, fontweight="bold")
    ax.text(shared_x, top_y, f"Shared aquatic foundation\n(n = {shared_n})", ha="center", va="bottom",
            color="#4D4D4D", fontsize=10, fontweight="bold")
    ax.text(aquatic_axis_x, top_y, f"Binary aquatic-dependence model\n(n = {aquatic_n})", ha="center", va="bottom",
            color="#006D5B", fontsize=9.3, fontweight="bold")

    # Local coefficient scale under each model axis.
    for axis_x, label in [(marine_axis_x, "Marine coefficient"), (aquatic_axis_x, "Aquatic coefficient")]:
        ax.hlines(-3.0, axis_x - 0.86, axis_x + 0.86, color="#333333", lw=0.8)
        ax.vlines(axis_x, -3.35, -2.65, color="#333333", lw=0.8)
        ax.text(axis_x - 0.86, -4.4, f"-{max_abs:.2f}", ha="center", va="top", fontsize=7)
        ax.text(axis_x, -4.4, "0", ha="center", va="top", fontsize=7)
        ax.text(axis_x + 0.86, -4.4, f"{max_abs:.2f}", ha="center", va="top", fontsize=7)
        ax.text(axis_x, -6.2, label, ha="center", va="top", fontsize=8.5)

    ax.text(-4.05, top_y - 1.2, "positive coefficient\nred, higher-rate", ha="left", va="top",
            fontsize=7.5, color=RED)
    ax.text(-4.05, top_y - 7.2, "negative coefficient\nblue, lower-rate", ha="left", va="top",
            fontsize=7.5, color=BLUE)

    # Legends. Keep them as basic matplotlib objects so the SVG remains editable.
    gene_set_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=SHARED_GREY, markeredgecolor="none",
               markersize=6, label="Shared aquatic foundation"),
        Line2D([0], [0], marker="^", color="none", markerfacecolor=MARINE_BLUE, markeredgecolor="none",
               markersize=6, label="Marine-specific predictors"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=AQUATIC_TEAL, markeredgecolor="none",
               markersize=6, label="Aquatic-specific predictors"),
    ]
    sign_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=RED, markeredgecolor="none",
               markersize=7, label="Positive coefficient"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=BLUE, markeredgecolor="none",
               markersize=7, label="Negative coefficient"),
    ]
    size_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", markeredgecolor="none",
               markersize=3.3, label="|coef| <= 0.1"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#555555", markeredgecolor="none",
               markersize=6.3, label="|coef| > 0.1"),
    ]
    leg1 = ax.legend(handles=gene_set_handles, title="Gene set", loc="upper left",
                     bbox_to_anchor=(1.01, 0.72), frameon=False, fontsize=8.5, title_fontsize=9.5)
    ax.add_artist(leg1)
    leg2 = ax.legend(handles=sign_handles, title="Coefficient sign", loc="upper left",
                     bbox_to_anchor=(1.01, 0.49), frameon=False, fontsize=8.5, title_fontsize=9.5)
    ax.add_artist(leg2)
    ax.legend(handles=size_handles, title="Coefficient size", loc="upper left",
              bbox_to_anchor=(1.01, 0.30), frameon=False, fontsize=8.5, title_fontsize=9.5)

    ax.set_title("Coefficient architecture of LASSO predictors", loc="left",
                 fontsize=13, fontweight="bold", pad=10)
    ax.set_xlim(-4.2, 4.25)
    ax.set_ylim(-7.5, top_y + 4)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.04, right=0.78, top=0.92, bottom=0.10)
    fig.savefig(PDF_FILE, transparent=True)
    fig.savefig(SVG_FILE, transparent=True)
    fig.savefig(PNG_FILE, transparent=True, dpi=300)
    plt.close(fig)

    print("Saved:")
    print(PDF_FILE)
    print(SVG_FILE)
    print(PNG_FILE)
    print(PLOT_TABLE_FILE)


if __name__ == "__main__":
    main()
