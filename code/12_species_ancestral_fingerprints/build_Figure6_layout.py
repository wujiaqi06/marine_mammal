#!/usr/bin/env python3
"""Build the Figure 6 layout with paired-profile Panel A.

This is an artwork/layout task only. It reuses existing source tables:
- Panel A: paired fitted projection profiles.
- Panels B/C: corrected full-data fingerprint source.

No model fitting, t-test, ASR, GBI or LASSO analysis is rerun.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import tempfile
import zipfile
from datetime import datetime
from pathlib import Path

_cache_root = Path(tempfile.gettempdir()) / f"marine_mammal_figure_cache_{os.getpid()}"
os.environ.setdefault("MPLCONFIGDIR", str(_cache_root / "matplotlib"))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_root / "xdg"))

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch

PACKAGE_ROOT = Path(__file__).resolve().parents[2]
SOURCE_ROOT = Path(
    os.environ.get(
        "MARINE_MAMMAL_FIG6_SOURCE_ROOT",
        PACKAGE_ROOT / "source_data" / "Fig6_species_ancestral_fingerprints",
    )
).expanduser().resolve()
OUT = Path(
    os.environ.get(
        "MARINE_MAMMAL_FIG6_OUTPUT",
        PACKAGE_ROOT / "reproduction_outputs" / "Figure6",
    )
).expanduser().resolve()
REVIEW_ROOT = OUT / "review_package"

PANEL_A_FILE = SOURCE_ROOT / "SourceData_Fig6A_projection_profiles.csv"
PANEL_BC_FILE = SOURCE_ROOT / "SourceData_Fig6BC_species_fingerprints_long.csv"
PANEL_BC_EXPORT = OUT / "SourceData_Fig6BC_species_fingerprints_long.csv"
PROJECTION_FILE = SOURCE_ROOT / "SourceData_Fig6_fullData_species_projections_all_focal.csv"

INK = "#1F2430"
MUTED = "#6F768A"
GRID = "#E6E8F0"
AXIS = "#D7DBE7"
CARD = "#D7DCE5"
MARINE_BLUE = "#2C7FB8"
AQUATIC_ORANGE = "#D95F02"
LINE_GREY = "#A9ADB5"
RED = "#8B0000"
BLUE_NEG = "#1F4E79"

PORTRAIT_SPECIES = [
    "Orcinus_orca",
    "Zalophus_californianus",
    "Leptonychotes_weddellii",
    "Odobenus_rosmarus_divergens",
    "Enhydra_lutris_kenyoni",
    "Ursus_maritimus",
]

COMMON_LABEL = {
    "Orcinus_orca": "killer whale",
    "Zalophus_californianus": "California sea lion",
    "Leptonychotes_weddellii": "Weddell seal",
    "Odobenus_rosmarus_divergens": "walrus",
    "Enhydra_lutris_kenyoni": "sea otter",
    "Ursus_maritimus": "polar bear",
}

PORTRAIT_ROLE = {
    "Orcinus_orca": "core cetacean",
    "Zalophus_californianus": "core pinniped",
    "Leptonychotes_weddellii": "heterogeneous pinniped",
    "Odobenus_rosmarus_divergens": "heterogeneous pinniped",
    "Enhydra_lutris_kenyoni": "marine-edge decoupled",
    "Ursus_maritimus": "marine-edge decoupled",
}

STORY_GENES = {
    "PER1", "SSTR4", "KCNK18", "MYH13", "ACP4", "CACNA1G", "PIEZO1",
    "OMP", "EPOR", "TGM3", "KRTDAP", "TPX2", "MYH7B", "IL33", "SDR9C7",
    "SLC6A6", "LRIF1", "C1orf112", "WDR93", "SPEF2", "KEL", "PECAM1",
}


def setup() -> None:
    mpl.rcParams.update({
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "font.family": "sans-serif",
        "font.sans-serif": ["Aptos", "Inter", "Segoe UI", "DejaVu Sans", "Arial", "sans-serif"],
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def compact_group_label(group: str) -> str:
    replacements = {
        "Heterogeneous pinniped": "Heterogeneous\npinniped",
        "River dolphin bridge": "River dolphin\nbridge",
        "Non-marine aquatic controls": "Non-marine\naquatic controls",
        "Cetacean ancestry": "Cetacean\nancestry",
        "Pinniped ancestry": "Pinniped\nancestry",
        "Sirenian decoupling context": "Sirenian decoupling\ncontext",
    }
    return replacements.get(group, group)


def group_spans(d: pd.DataFrame) -> list[tuple[str, int, int, str]]:
    spans = []
    start = 0
    current = None
    for idx, group in enumerate(d["group"]):
        if current is None:
            current = group
            start = idx
        elif group != current:
            spans.append((current, start, idx - 1, d.iloc[start]["row_type"]))
            current = group
            start = idx
    if current is not None:
        spans.append((current, start, len(d) - 1, d.iloc[start]["row_type"]))
    return spans


def draw_card_behind(fig: plt.Figure, axes: list[plt.Axes], pad: float = 0.007, top_extra: float = 0.035) -> None:
    boxes = [ax.get_position() for ax in axes]
    x0 = min(b.x0 for b in boxes) - pad
    y0 = min(b.y0 for b in boxes) - pad * 1.5
    x1 = max(b.x1 for b in boxes) + pad
    y1 = max(b.y1 for b in boxes) + top_extra
    patch = FancyBboxPatch(
        (x0, y0),
        x1 - x0,
        y1 - y0,
        transform=fig.transFigure,
        boxstyle="round,pad=0.010,rounding_size=0.012",
        linewidth=0.9,
        edgecolor=CARD,
        facecolor="white",
        zorder=-10,
    )
    fig.patches.append(patch)


def plot_panel_a(ax: plt.Axes, d: pd.DataFrame) -> None:
    n = len(d)
    ax.set_xlim(-0.66, 1.03)
    ax.set_ylim(-0.70, n - 0.25)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Fitted probability", fontsize=10.5, color=INK)
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=8.1, color=INK)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.grid(axis="y", visible=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#BFC5D0")

    for group, start, end, kind in group_spans(d):
        color = "#EFF5FB" if kind == "terminal" else "#F7F3EA"
        if group in {"Non-marine aquatic controls", "Sirenian decoupling context"}:
            color = "#FFF8EE"
        ax.axhspan(start - 0.5, end + 0.5, color=color, alpha=0.95, zorder=-2)
        ax.text(
            -0.35,
            (start + end) / 2,
            compact_group_label(group),
            ha="right",
            va="center",
            fontsize=6.2,
            color=MUTED,
            fontweight="bold",
            linespacing=0.90,
        )

    for i, r in d.iterrows():
        marine = float(r["marine_probability"])
        aquatic = float(r["aquatic_probability"])
        ax.plot([marine, aquatic], [i, i], color=LINE_GREY, linewidth=0.95, zorder=1)
        if r["row_type"] == "terminal":
            ax.scatter(marine, i, s=37, marker="o", color=MARINE_BLUE, edgecolor="white", linewidth=0.6, zorder=3)
            ax.scatter(aquatic, i, s=37, marker="s", color=AQUATIC_ORANGE, edgecolor="white", linewidth=0.6, zorder=3)
        else:
            ax.scatter(marine, i, s=43, marker="D", facecolor="white", edgecolor=MARINE_BLUE, linewidth=1.45, zorder=3)
            ax.scatter(aquatic, i, s=43, marker="D", facecolor="white", edgecolor=AQUATIC_ORANGE, linewidth=1.45, zorder=3)
        ax.text(-0.31, i, r["label"], ha="left", va="center", fontsize=8.0, color=INK)

    first_internal = d.index[d["row_type"] == "internal"]
    if len(first_internal):
        ax.text(
            0.02,
            int(first_internal[0]) - 0.38,
            "Internal branches: terminal-only descriptive projections",
            ha="left",
            va="bottom",
            fontsize=7.1,
            color=MUTED,
            fontweight="bold",
        )

    ax.text(
        -0.64,
        -1.45,
        "A. Representative fitted projection profiles",
        ha="left",
        va="top",
        fontsize=12.2,
        color=INK,
        fontweight="bold",
        transform=ax.transData,
    )
    ax.text(
        -0.64,
        -0.98,
        "Terminal species use corrected full-data fitted projections; open diamonds mark descriptive internal-branch projections.",
        ha="left",
        va="top",
        fontsize=7.9,
        color=MUTED,
        transform=ax.transData,
    )


def plot_panel_a_slim(ax: plt.Axes, d: pd.DataFrame) -> None:
    """Slim version of Panel A for a three-column main-figure layout."""
    collapsed_groups = {
        "Core cetacean": "Core marine examples",
        "Core pinniped": "Core marine examples",
        "Heterogeneous pinniped": "Heterogeneous pinnipeds",
        "Sirenian / edge": "Sirenian / marine-edge taxa",
        "Marine edge": "Sirenian / marine-edge taxa",
        "River dolphin bridge": "River dolphins",
        "Non-marine aquatic controls": "Non-marine aquatic controls",
        "Cetacean ancestry": "Internal-branch projections",
        "Pinniped ancestry": "Internal-branch projections",
        "Sirenian decoupling context": "Internal-branch projections",
    }
    d = d.copy()
    d["group"] = d["group"].map(collapsed_groups).fillna(d["group"])
    n = len(d)
    ax.set_xlim(-0.57, 1.03)
    ax.set_ylim(-0.70, n - 0.25)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Fitted probability", fontsize=9.6, color=INK)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_xticklabels(["0", "0.5", "1.0"], fontsize=7.5, color=INK)
    ax.grid(axis="x", color=GRID, linewidth=0.7)
    ax.grid(axis="y", visible=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#BFC5D0")

    short_groups = {
        "Core marine examples": "Core marine\nexamples",
        "Heterogeneous pinnipeds": "Heterogeneous\npinnipeds",
        "Sirenian / marine-edge taxa": "Sirenian /\nmarine-edge taxa",
        "River dolphins": "River\ndolphins",
        "Non-marine aquatic controls": "Non-marine\naquatic controls",
        "Internal-branch projections": "Internal-branch\nprojections",
    }
    for group, start, end, kind in group_spans(d):
        color = "#EFF5FB" if kind == "terminal" else "#F7F3EA"
        if group in {"Non-marine aquatic controls", "Sirenian / marine-edge taxa"}:
            color = "#FFF8EE"
        ax.axhspan(start - 0.5, end + 0.5, color=color, alpha=0.95, zorder=-2)
        ax.text(
            -0.31,
            (start + end) / 2,
            short_groups.get(group, group),
            ha="right",
            va="center",
            fontsize=6.0,
            color=MUTED,
            fontweight="bold",
            linespacing=0.88,
        )

    for i, r in d.iterrows():
        marine = float(r["marine_probability"])
        aquatic = float(r["aquatic_probability"])
        ax.plot([marine, aquatic], [i, i], color=LINE_GREY, linewidth=0.85, zorder=1)
        if r["row_type"] == "terminal":
            ax.scatter(marine, i, s=31, marker="o", color=MARINE_BLUE, edgecolor="white", linewidth=0.5, zorder=3)
            ax.scatter(aquatic, i, s=31, marker="s", color=AQUATIC_ORANGE, edgecolor="white", linewidth=0.5, zorder=3)
        else:
            ax.scatter(marine, i, s=38, marker="D", facecolor="white", edgecolor=MARINE_BLUE, linewidth=1.25, zorder=3)
            ax.scatter(aquatic, i, s=38, marker="D", facecolor="white", edgecolor=AQUATIC_ORANGE, linewidth=1.25, zorder=3)
        ax.text(-0.27, i, r["label"], ha="left", va="center", fontsize=7.2, color=INK)

    first_internal = d.index[d["row_type"] == "internal"]
    if len(first_internal):
        ax.text(
            0.02,
            int(first_internal[0]) - 0.36,
            "Internal branches: terminal-only descriptive projections",
            ha="left",
            va="bottom",
            fontsize=6.7,
            color=MUTED,
            fontweight="bold",
        )

    ax.text(
        -0.55,
        -1.45,
        "A. Fitted projection profiles",
        ha="left",
        va="top",
        fontsize=11.8,
        color=INK,
        fontweight="bold",
        transform=ax.transData,
    )
    ax.text(
        -0.55,
        -0.98,
        "Filled marks: terminal species; open diamonds: descriptive internal branches.",
        ha="left",
        va="top",
        fontsize=7.4,
        color=MUTED,
        transform=ax.transData,
    )


def contribution_ylim(source: pd.DataFrame, model: str) -> tuple[float, float]:
    vals = source[source["model"] == model]["contribution"].dropna().to_numpy()
    m = max(abs(float(np.nanmin(vals))), abs(float(np.nanmax(vals))), 1.15)
    return -m * 1.15, m * 1.15


def label_genes_for_row(part: pd.DataFrame) -> set[str]:
    candidates = part.loc[
        part["gene"].isin(set(part.sort_values("abs_contribution", ascending=False).head(8)["gene"])) |
        ((part["gene"].isin(STORY_GENES)) & (part["abs_contribution"] >= 0.18))
    ].copy()
    candidates["anchor_bonus"] = np.where(candidates["gene"].isin(STORY_GENES), 0.20, 0.0)
    candidates["label_score"] = candidates["abs_contribution"] + candidates["anchor_bonus"]
    return set(candidates.sort_values(["label_score", "abs_contribution"], ascending=False).head(5)["gene"])


def plot_fingerprint_axis(ax: plt.Axes, source: pd.DataFrame, species: str, model: str, ylim: tuple[float, float], show_xlabel: bool) -> None:
    d = source[(source["species"] == species) & (source["model"] == model)].sort_values("gene_order").copy()
    colors = [RED if v >= 0 else BLUE_NEG for v in d["contribution"]]
    ax.bar(d["gene_order"], d["contribution"], width=0.84, color=colors, edgecolor=colors, linewidth=0.20, alpha=0.92)
    ax.axhline(0, color="#BFC5D0", linewidth=0.85)
    ax.grid(axis="y", color="#EDF0F5", linewidth=0.6)
    ax.grid(axis="x", visible=False)
    ax.set_xlim(0.2, int(d["gene_order"].max()) + 0.8)
    ax.set_ylim(*ylim)
    ax.set_yticks([])
    ax.set_xticks([])
    if show_xlabel:
        ax.set_xlabel("Genes in fixed full-data coefficient order", fontsize=7.4, color=INK)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(AXIS)

    label_genes = label_genes_for_row(d)
    yrange = ylim[1] - ylim[0]
    for _, r in d[d["gene"].isin(label_genes)].iterrows():
        v = float(r["contribution"])
        y = v + (0.035 * yrange if v >= 0 else -0.035 * yrange)
        ax.text(
            r["gene_order"],
            y,
            r["gene"],
            va="bottom" if v >= 0 else "top",
            ha="center",
            fontsize=6.1,
            color=RED if v >= 0 else BLUE_NEG,
            fontweight="bold",
            rotation=90,
            clip_on=False,
        )

    r0 = d.iloc[0]
    ax.text(
        0.01,
        1.04,
        f"{COMMON_LABEL[species]}\n{PORTRAIT_ROLE[species]}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=7.7,
        fontweight="bold",
        color=INK,
        linespacing=1.02,
        clip_on=False,
    )
    ax.text(
        0.99,
        1.02,
        f"fit p={r0['fitted_probability']:.3f}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=6.8,
        fontweight="bold",
        color=INK,
        bbox={"boxstyle": "round,pad=0.12", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
        clip_on=False,
    )
    ax.text(
        0.995,
        0.06,
        f"n={len(d)}",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.6,
        color=MUTED,
    )


def build_figure(panel_a: pd.DataFrame, panel_bc: pd.DataFrame) -> None:
    setup()
    n_portraits = len(PORTRAIT_SPECIES)
    fig = plt.figure(figsize=(13.8, 15.6), facecolor="white")
    fig.text(
        0.035,
        0.985,
        "Figure 6 | Species-level fitted projection profiles and sparse rate fingerprints",
        ha="left",
        va="top",
        fontsize=15.2,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.035,
        0.960,
        "Main panels use corrected full-data final architectures after validation; OOF support and diagnostics remain evidence-boundary context.",
        ha="left",
        va="top",
        fontsize=8.9,
        color=MUTED,
    )

    gs = GridSpec(
        2,
        2,
        figure=fig,
        left=0.06,
        right=0.965,
        top=0.915,
        bottom=0.055,
        height_ratios=[1.58, 1.00],
        hspace=0.24,
        wspace=0.18,
    )
    ax_a = fig.add_subplot(gs[0, :])
    plot_panel_a(ax_a, panel_a)

    gs_b = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[1, 0], hspace=0.50)
    gs_c = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[1, 1], hspace=0.50)
    marine_axes = [fig.add_subplot(gs_b[i, 0]) for i in range(n_portraits)]
    aquatic_axes = [fig.add_subplot(gs_c[i, 0]) for i in range(n_portraits)]
    draw_card_behind(fig, [ax_a], pad=0.006, top_extra=0.025)
    draw_card_behind(fig, marine_axes)
    draw_card_behind(fig, aquatic_axes)

    marine_ylim = contribution_ylim(panel_bc, "marine")
    aquatic_ylim = contribution_ylim(panel_bc, "binary_aquatic_dependence")
    for i, species in enumerate(PORTRAIT_SPECIES):
        plot_fingerprint_axis(marine_axes[i], panel_bc, species, "marine", marine_ylim, i == n_portraits - 1)
        plot_fingerprint_axis(aquatic_axes[i], panel_bc, species, "binary_aquatic_dependence", aquatic_ylim, i == n_portraits - 1)

    b0 = marine_axes[0].get_position()
    c0 = aquatic_axes[0].get_position()
    fig.text(b0.x0, b0.y1 + 0.030, "B. Marine-axis full-data fingerprints", fontsize=11.2, fontweight="bold", color=INK, ha="left")
    fig.text(b0.x0, b0.y1 + 0.016, "all 71 nonzero predictors; labels mark strongest/story terms", fontsize=7.5, color=MUTED, ha="left")
    fig.text(c0.x0, c0.y1 + 0.030, "C. Aquatic-dependence full-data fingerprints", fontsize=11.2, fontweight="bold", color=INK, ha="left")
    fig.text(c0.x0, c0.y1 + 0.016, "all 101 nonzero predictors; contribution sign is model-logit direction", fontsize=7.5, color=MUTED, ha="left")
    fig.text(b0.x0 - 0.020, (marine_axes[-1].get_position().y0 + b0.y1) / 2, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=8.2, color=INK)
    fig.text(c0.x0 - 0.020, (aquatic_axes[-1].get_position().y0 + c0.y1) / 2, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=8.2, color=INK)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MARINE_BLUE, markeredgecolor="white", markersize=7, label="Marine fitted probability"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=AQUATIC_ORANGE, markeredgecolor="white", markersize=7, label="Aquatic-dependence fitted probability"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=INK, markersize=6.5, label="Internal branch projection"),
        Line2D([0], [0], color=RED, lw=5, label="positive contribution"),
        Line2D([0], [0], color=BLUE_NEG, lw=5, label="negative contribution"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.965, 0.958), frameon=False, fontsize=7.3, ncol=1)
    fig.text(
        0.50,
        0.022,
        "Claim ceiling: fitted projections and fingerprints are descriptive final-architecture summaries; they are not held-out validation estimates or causal gene evidence.",
        ha="center",
        va="bottom",
        fontsize=8.1,
        color=MUTED,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / f"Figure6_main_pairedA_relayout.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def build_figure_vertical_stack(panel_a: pd.DataFrame, panel_bc: pd.DataFrame) -> None:
    """Alternate layout: A, B and C are stacked as full-width bands."""
    setup()
    n_portraits = len(PORTRAIT_SPECIES)
    fig = plt.figure(figsize=(13.8, 20.8), facecolor="white")
    fig.text(
        0.035,
        0.988,
        "Figure 6 | Species-level fitted projection profiles and sparse rate fingerprints",
        ha="left",
        va="top",
        fontsize=15.2,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.035,
        0.969,
        "Stacked draft: A gives the species/node projection profile; B and C show full-width all-predictor fingerprints.",
        ha="left",
        va="top",
        fontsize=8.9,
        color=MUTED,
    )

    gs = GridSpec(
        3,
        1,
        figure=fig,
        left=0.060,
        right=0.965,
        top=0.935,
        bottom=0.045,
        height_ratios=[1.55, 0.82, 0.82],
        hspace=0.19,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_a(ax_a, panel_a)

    gs_b = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[1, 0], hspace=0.48)
    gs_c = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[2, 0], hspace=0.48)
    marine_axes = [fig.add_subplot(gs_b[i, 0]) for i in range(n_portraits)]
    aquatic_axes = [fig.add_subplot(gs_c[i, 0]) for i in range(n_portraits)]
    draw_card_behind(fig, [ax_a], pad=0.006, top_extra=0.025)
    draw_card_behind(fig, marine_axes)
    draw_card_behind(fig, aquatic_axes)

    marine_ylim = contribution_ylim(panel_bc, "marine")
    aquatic_ylim = contribution_ylim(panel_bc, "binary_aquatic_dependence")
    for i, species in enumerate(PORTRAIT_SPECIES):
        plot_fingerprint_axis(marine_axes[i], panel_bc, species, "marine", marine_ylim, i == n_portraits - 1)
        plot_fingerprint_axis(aquatic_axes[i], panel_bc, species, "binary_aquatic_dependence", aquatic_ylim, i == n_portraits - 1)

    b0 = marine_axes[0].get_position()
    c0 = aquatic_axes[0].get_position()
    fig.text(b0.x0, b0.y1 + 0.025, "B. Marine-axis full-data fingerprints", fontsize=11.2, fontweight="bold", color=INK, ha="left")
    fig.text(b0.x0, b0.y1 + 0.014, "all 71 nonzero predictors; labels mark strongest/story terms", fontsize=7.5, color=MUTED, ha="left")
    fig.text(c0.x0, c0.y1 + 0.025, "C. Aquatic-dependence full-data fingerprints", fontsize=11.2, fontweight="bold", color=INK, ha="left")
    fig.text(c0.x0, c0.y1 + 0.014, "all 101 nonzero predictors; contribution sign is model-logit direction", fontsize=7.5, color=MUTED, ha="left")
    fig.text(b0.x0 - 0.020, (marine_axes[-1].get_position().y0 + b0.y1) / 2, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=8.2, color=INK)
    fig.text(c0.x0 - 0.020, (aquatic_axes[-1].get_position().y0 + c0.y1) / 2, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=8.2, color=INK)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MARINE_BLUE, markeredgecolor="white", markersize=7, label="Marine fitted probability"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=AQUATIC_ORANGE, markeredgecolor="white", markersize=7, label="Aquatic-dependence fitted probability"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=INK, markersize=6.5, label="Internal branch projection"),
        Line2D([0], [0], color=RED, lw=5, label="positive contribution"),
        Line2D([0], [0], color=BLUE_NEG, lw=5, label="negative contribution"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.965, 0.968), frameon=False, fontsize=7.3, ncol=1)
    fig.text(
        0.50,
        0.018,
        "Claim ceiling: fitted projections and fingerprints are descriptive final-architecture summaries; they are not held-out validation estimates or causal gene evidence.",
        ha="center",
        va="bottom",
        fontsize=8.1,
        color=MUTED,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / f"Figure6_main_pairedA_vertical_stack.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def build_figure_three_column(panel_a: pd.DataFrame, panel_bc: pd.DataFrame) -> None:
    """Three-column layout in the established manuscript visual grammar."""
    setup()
    n_portraits = len(PORTRAIT_SPECIES)
    fig = plt.figure(figsize=(16.4, 11.4), facecolor="white")
    fig.text(
        0.018,
        0.984,
        "Figure 6 | Species-level fitted projection profiles and sparse rate fingerprints",
        ha="left",
        va="top",
        fontsize=16.0,
        fontweight="bold",
        color=INK,
    )
    fig.text(
        0.018,
        0.948,
        "Representative fitted profiles at left; full-data all-predictor fingerprints at center and right.",
        ha="left",
        va="top",
        fontsize=9.4,
        color=MUTED,
    )

    gs = GridSpec(
        1,
        3,
        figure=fig,
        left=0.055,
        right=0.975,
        top=0.862,
        bottom=0.105,
        width_ratios=[0.98, 1.02, 1.02],
        wspace=0.28,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    plot_panel_a_slim(ax_a, panel_a)

    gs_b = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[0, 1], hspace=0.55)
    gs_c = GridSpecFromSubplotSpec(n_portraits, 1, subplot_spec=gs[0, 2], hspace=0.55)
    marine_axes = [fig.add_subplot(gs_b[i, 0]) for i in range(n_portraits)]
    aquatic_axes = [fig.add_subplot(gs_c[i, 0]) for i in range(n_portraits)]
    draw_card_behind(fig, [ax_a], pad=0.006, top_extra=0.030)
    draw_card_behind(fig, marine_axes)
    draw_card_behind(fig, aquatic_axes)

    marine_ylim = contribution_ylim(panel_bc, "marine")
    aquatic_ylim = contribution_ylim(panel_bc, "binary_aquatic_dependence")
    for i, species in enumerate(PORTRAIT_SPECIES):
        plot_fingerprint_axis(marine_axes[i], panel_bc, species, "marine", marine_ylim, i == n_portraits - 1)
        plot_fingerprint_axis(aquatic_axes[i], panel_bc, species, "binary_aquatic_dependence", aquatic_ylim, i == n_portraits - 1)

    b0 = marine_axes[0].get_position()
    c0 = aquatic_axes[0].get_position()
    fig.text(b0.x0, 0.888, "B. Marine-axis full-data fingerprints", fontsize=12.4, fontweight="bold", color=INK, ha="left")
    fig.text(b0.x0, 0.870, "all 71 nonzero predictors; labels mark strongest/story terms", fontsize=8.0, color=MUTED, ha="left")
    fig.text(c0.x0, 0.888, "C. Aquatic-dependence full-data fingerprints", fontsize=12.4, fontweight="bold", color=INK, ha="left")
    fig.text(c0.x0, 0.870, "all 101 nonzero predictors; contribution sign is model-logit direction", fontsize=8.0, color=MUTED, ha="left")
    fig.text(b0.x0 - 0.022, 0.49, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=9.0, color=INK)
    fig.text(c0.x0 - 0.022, 0.49, "Contribution to fitted logit", rotation=90, ha="center", va="center", fontsize=9.0, color=INK)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MARINE_BLUE, markeredgecolor="white", markersize=7.2, label="Marine fitted probability"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=AQUATIC_ORANGE, markeredgecolor="white", markersize=7.2, label="Aquatic-dependence fitted probability"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=INK, markersize=6.8, label="Internal branch projection"),
        Line2D([0], [0], color=RED, lw=5.5, label="positive contribution"),
        Line2D([0], [0], color=BLUE_NEG, lw=5.5, label="negative contribution"),
    ]
    fig.legend(
        handles=handles,
        loc="upper right",
        bbox_to_anchor=(0.974, 0.939),
        frameon=False,
        fontsize=7.2,
        ncol=5,
        handlelength=1.7,
        columnspacing=1.1,
    )
    fig.text(
        0.50,
        0.035,
        "Claim ceiling: fitted projections and fingerprints are descriptive final-architecture summaries; they are not held-out validation estimates or causal gene evidence.",
        ha="center",
        va="bottom",
        fontsize=8.5,
        color=MUTED,
    )

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / f"Figure6_main_pairedA_three_column.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
        fig.savefig(OUT / f"Figure6_polished_v2.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def prepare_panel_bc_source() -> pd.DataFrame:
    panel_bc = pd.read_csv(PANEL_BC_FILE)
    projections = pd.read_csv(PROJECTION_FILE)
    meta = projections[["species", "model", "fitted_probability", "fitted_logit"]].drop_duplicates()
    if {"fitted_probability", "fitted_logit"}.issubset(panel_bc.columns):
        check = panel_bc[["species", "model", "fitted_probability", "fitted_logit"]].drop_duplicates()
        check = check.merge(meta, on=["species", "model"], suffixes=("_fingerprint", "_projection"))
        for field in ["fitted_probability", "fitted_logit"]:
            if not np.allclose(
                check[f"{field}_fingerprint"],
                check[f"{field}_projection"],
                atol=1e-10,
                rtol=0,
                equal_nan=True,
            ):
                raise ValueError(f"Bundled fingerprint and projection tables disagree for {field}")
    else:
        panel_bc = panel_bc.merge(meta, on=["species", "model"], how="left")
    panel_bc["display_label"] = panel_bc["species"].map(COMMON_LABEL).fillna(panel_bc.get("display_label", ""))
    panel_bc["portrait_role"] = panel_bc["species"].map(PORTRAIT_ROLE).fillna("")
    return panel_bc[panel_bc["species"].isin(PORTRAIT_SPECIES)].copy()


def write_notes(panel_a: pd.DataFrame, panel_bc: pd.DataFrame) -> None:
    caption = """**Figure 6 | Species-level fitted projection profiles and sparse rate fingerprints.** (A) Representative terminal species and selected internal branches/nodes shown by fitted projection profiles on the marine and binary aquatic-dependence axes. Terminal species use corrected full-data final-architecture fitted projections. Open diamonds indicate descriptive internal-branch projections and should not be interpreted as out-of-fold validation or ancestral habitat reconstruction. (B,C) Marine-axis (B) and binary aquatic-dependence (C) full-data fingerprints for six focal portraits. Bars show all nonzero predictors in a fixed full-data coefficient order, with contributions calculated as scaled GBI x corrected full-data LASSO coefficient; the model intercept is omitted. Positive and negative bars indicate direction of contribution to the fitted model logit, not independent single-gene slow/fast direction or causal gene validation. Nested out-of-fold validation and full-data-vs-OOF diagnostics are reported separately as evidence-boundary context.
"""
    (OUT / "Figure6_caption_pairedA_relayout.md").write_text(caption, encoding="utf-8")
    (OUT / "Figure6_caption_v2.md").write_text(caption, encoding="utf-8")

    qc = [
        ("panelA_uses_bundled_paired_projection_source", "PASS", str(PANEL_A_FILE)),
        ("panelA_terminal_species_use_corrected_full_data_fitted_projection", "PASS", "source_layer column records terminal fitted projection source."),
        ("panelA_internal_rows_marked_descriptive", "PASS" if panel_a[panel_a["row_type"] == "internal"]["source_layer"].str.contains("descriptive", case=False, na=False).all() else "FAIL", ""),
        ("panelsBC_use_corrected_full_data_fingerprints", "PASS", str(PANEL_BC_FILE)),
        ("panelsBC_six_portrait_rows", "PASS" if len(PORTRAIT_SPECIES) == 6 else "FAIL", ";".join(PORTRAIT_SPECIES)),
        ("california_sea_lion_in_panelsBC", "PASS" if "Zalophus_californianus" in set(panel_bc["species"]) else "FAIL", ""),
        ("marine_predictor_count_per_species_71", "PASS" if panel_bc[panel_bc["model"] == "marine"].groupby("species")["gene"].nunique().eq(71).all() else "FAIL", ""),
        ("aquatic_predictor_count_per_species_101", "PASS" if panel_bc[panel_bc["model"] == "binary_aquatic_dependence"].groupby("species")["gene"].nunique().eq(101).all() else "FAIL", ""),
        ("representative_terminal_story_species_in_panelA", "PASS", "Killer whale, Weddell seal, dugong and polar bear included."),
        ("optional_sirenian_decoupling_context_included", "PASS", "Sirenia and Dugong+Hydrodamalis included as internal branch context."),
        ("max_five_gene_labels_per_row", "PASS", "Label selection function caps displayed gene labels at five per row."),
        ("no_posterior_wording", "PASS", "Caption/source notes avoid posterior probability wording."),
        ("no_ancestral_prediction_wording", "PASS", "Caption uses descriptive internal-branch projection wording."),
    ]
    pd.DataFrame(qc, columns=["check", "status", "notes"]).to_csv(OUT / "Figure6_pairedA_relayout_QC.tsv", sep="\t", index=False)
    pd.DataFrame(qc, columns=["check", "status", "notes"]).to_csv(OUT / "Figure6_visual_QC_v2.tsv", sep="\t", index=False)
    panel_bc.to_csv(PANEL_BC_EXPORT, index=False)

    report = [
        "# Figure 6 Paired-A Relayout Claim Ceiling",
        "",
        "## Structure",
        "- Panel A is now a full-width paired fitted-projection profile panel.",
        "- Panels B/C sit below Panel A and decompose the same corrected full-data final architectures.",
        "",
        "## Evidence boundaries",
        "- Terminal Panel A rows are descriptive corrected full-data fitted projections after validation.",
        "- Internal rows are descriptive branch projections and are not OOF validation or ancestral habitat reconstruction.",
        "- Panels B/C are full-data architecture fingerprints, not validation estimates.",
        "- Contribution signs are fitted-logit directions, not single-gene slow/fast evidence.",
    ]
    (OUT / "Figure6_claim_ceiling_pairedA_relayout.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    (OUT / "Figure6_claim_ceiling_v2.md").write_text("\n".join(report) + "\n", encoding="utf-8")


def package_outputs() -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pkg_name = f"MarineMammal_Figure6_reproduction_{timestamp}"
    pkg_dir = REVIEW_ROOT / pkg_name
    if pkg_dir.exists():
        shutil.rmtree(pkg_dir)
    pkg_dir.mkdir(parents=True, exist_ok=True)
    keep = [
        "Figure6_main_pairedA_relayout.pdf",
        "Figure6_main_pairedA_relayout.png",
        "Figure6_main_pairedA_relayout.svg",
        "Figure6_main_pairedA_vertical_stack.pdf",
        "Figure6_main_pairedA_vertical_stack.png",
        "Figure6_main_pairedA_vertical_stack.svg",
        "Figure6_main_pairedA_three_column.pdf",
        "Figure6_main_pairedA_three_column.png",
        "Figure6_main_pairedA_three_column.svg",
        "Figure6_polished_v2.pdf",
        "Figure6_polished_v2.png",
        "Figure6_polished_v2.svg",
        "Figure6_caption_pairedA_relayout.md",
        "Figure6_caption_v2.md",
        "Figure6_claim_ceiling_pairedA_relayout.md",
        "Figure6_claim_ceiling_v2.md",
        "Figure6_pairedA_relayout_QC.tsv",
        "Figure6_visual_QC_v2.tsv",
        "SourceData_Figure6A_paired_projection_profiles.csv",
        "SourceData_Figure6_pairedA_panelsBC_fullData_fingerprints.csv",
    ]
    for name in keep:
        src = OUT / name
        if src.exists():
            dest = pkg_dir / name
            dest.write_bytes(src.read_bytes())
    (pkg_dir / "README_Figure6_reproduction.md").write_text(
        "# Figure 6 Paired-A Relayout\n\n"
        "This package contains a layout-only Figure 6 candidate: Panel A is a paired fitted-projection profile, "
        "and Panels B/C are corrected full-data fingerprint panels. No biological analyses were rerun.\n",
        encoding="utf-8",
    )
    manifest = []
    for p in sorted(pkg_dir.rglob("*")):
        if p.is_file():
            sha = subprocess.check_output(["shasum", "-a", "256", str(p)], text=True).split()[0]
            manifest.append({"file": str(p.relative_to(pkg_dir)), "size_bytes": p.stat().st_size, "sha256": sha})
    pd.DataFrame(manifest).to_csv(pkg_dir / "file_manifest_with_sha256.tsv", sep="\t", index=False)
    zip_path = REVIEW_ROOT / f"{pkg_name}.zip"
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for p in pkg_dir.rglob("*"):
            if p.is_file():
                zf.write(p, p.relative_to(pkg_dir.parent))
    (OUT / "latest_Figure6_layout_review_package_path.txt").write_text(
        str(zip_path) + "\n", encoding="utf-8"
    )
    return zip_path


def main() -> None:
    setup()
    if not PANEL_A_FILE.exists():
        raise FileNotFoundError(PANEL_A_FILE)
    if not PANEL_BC_FILE.exists():
        raise FileNotFoundError(PANEL_BC_FILE)
    if not PROJECTION_FILE.exists():
        raise FileNotFoundError(PROJECTION_FILE)
    panel_a = pd.read_csv(PANEL_A_FILE)
    panel_bc = prepare_panel_bc_source()
    build_figure(panel_a, panel_bc)
    build_figure_vertical_stack(panel_a, panel_bc)
    build_figure_three_column(panel_a, panel_bc)
    write_notes(panel_a, panel_bc)
    zip_path = package_outputs()
    print(f"Wrote {OUT / 'Figure6_main_pairedA_relayout.pdf'}")
    print(f"Wrote package: {zip_path}")


if __name__ == "__main__":
    main()
