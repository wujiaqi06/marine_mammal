#!/usr/bin/env python3
"""
Figure 4C: Module-level partition of LASSO predictors.

Creates two matplotlib variants:
1. Count-only bubble matrix.
2. Full gene-name bubble matrix with every module-annotated gene listed per non-empty cell.

This panel is a visual architecture summary. Full gene lists remain in Table 1
and Supplementary Table S5.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, HPacker, TextArea

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
INPUT_DIR = FIGURE_DIR / "inputs"
OUTPUT_DIR = FIGURE_DIR / "Figure4C_outputs"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MODULE_FILE = INPUT_DIR / "Table1_module_summary_working.tsv"
MAIN_CANDIDATE_FILE = INPUT_DIR / "Table1_main_text_candidate.tsv"
GENE_LEVEL_FILE = INPUT_DIR / "Table1_gene_level_working_table.tsv"
PLOT_TABLE_FILE = OUTPUT_DIR / "Figure4C_plot_table.tsv"

COUNT_ONLY_PDF = OUTPUT_DIR / "Figure4C_module_partition_count_only.pdf"
COUNT_ONLY_SVG = OUTPUT_DIR / "Figure4C_module_partition_count_only.svg"
COUNT_ONLY_PNG = OUTPUT_DIR / "Figure4C_module_partition_count_only.png"
EXAMPLES_PDF = OUTPUT_DIR / "Figure4C_module_partition_examples.pdf"
EXAMPLES_SVG = OUTPUT_DIR / "Figure4C_module_partition_examples.svg"
EXAMPLES_PNG = OUTPUT_DIR / "Figure4C_module_partition_examples.png"
FINAL_PDF = OUTPUT_DIR / "Figure4C_module_annotated_partition.pdf"
FINAL_SVG = OUTPUT_DIR / "Figure4C_module_annotated_partition.svg"
FINAL_PNG = OUTPUT_DIR / "Figure4C_module_annotated_partition.png"

SHARED_GREY = "#9E9E9E"
MARINE_BLUE = "#4C78A8"
AQUATIC_TEAL = "#44A08D"
FAST_RED = "#8B0000"
SLOW_BLUE = "#1F4E79"
DARK = "#111111"
GRID = "#E6E6E6"

MODULE_ORDER = [
    "epithelial barrier / keratinization / body-surface interface",
    "hematological / platelet / circulatory regulation",
    "immune / cytokine / inflammatory regulation",
    "sensory remodeling / olfactory or visual systems",
    "ion channels / mechanosensation / membrane signaling",
    "sperm / cilia / CatSper / reproductive-cell function",
    "DNA repair / chromatin / genome maintenance",
    "Skeletal, muscle and extracellular-matrix remodeling",
    "metabolism / mitochondrial / redox regulation",
    "fatty-acid / lipid metabolism",
    "endocrine / reproductive hormone systems",
    "transcriptional / systemic regulatory genes",
]

MODULE_LABELS = {
    "epithelial barrier / keratinization / body-surface interface":
        "Epithelial barrier / keratinization /\nbody-surface interface",
    "hematological / platelet / circulatory regulation":
        "Hematological / platelet /\ncirculatory regulation",
    "immune / cytokine / inflammatory regulation":
        "Immune / cytokine /\ninflammatory regulation",
    "sensory remodeling / olfactory or visual systems":
        "Sensory remodeling /\nolfactory or visual systems",
    "ion channels / mechanosensation / membrane signaling":
        "Ion channels / mechanosensation /\nmembrane signaling",
    "sperm / cilia / CatSper / reproductive-cell function":
        "Sperm / cilia / CatSper /\nreproductive-cell function",
    "DNA repair / chromatin / genome maintenance":
        "DNA repair / chromatin /\ngenome maintenance",
    "Skeletal, muscle and extracellular-matrix remodeling":
        "Skeletal, muscle and\nextracellular-matrix remodeling",
    "metabolism / mitochondrial / redox regulation":
        "Metabolism / mitochondrial /\nredox regulation",
    "fatty-acid / lipid metabolism":
        "Fatty-acid / lipid\nmetabolism",
    "endocrine / reproductive hormone systems":
        "Endocrine / reproductive\nhormone systems",
    "transcriptional / systemic regulatory genes":
        "Transcriptional / systemic\nregulatory genes",
}

CATEGORIES = [
    ("Shared aquatic foundation", "n_shared", "shared_genes", SHARED_GREY, 0.0),
    ("Marine-specific predictors", "n_marine_specific", "marine_specific_genes", MARINE_BLUE, 1.28),
    ("Aquatic-specific predictors", "n_aquatic_specific", "aquatic_specific_genes", AQUATIC_TEAL, 2.65),
]


def normalize_module_name(name: str) -> str:
    return name.replace("reproductive cell", "reproductive-cell")


def require_columns(df: pd.DataFrame, cols: list[str], label: str) -> None:
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in {label}: {', '.join(missing)}")


def full_gene_label(text: object) -> str:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return ""
    parts = [p.strip() for p in str(text).replace(";", ",").split(",") if p.strip()]
    if not parts:
        return ""
    return ", ".join(parts)


def split_genes(text: object) -> list[str]:
    if text is None or (isinstance(text, float) and np.isnan(text)):
        return []
    return [p.strip() for p in str(text).replace(";", ",").split(",") if p.strip()]


def coef_color(coef: float) -> str:
    return FAST_RED if coef >= 0 else SLOW_BLUE


def format_gene_label(gene: str, context: str, coef_lookup: dict[str, dict[str, float]]) -> dict[str, object]:
    vals = coef_lookup.get(gene)
    if vals is None:
        return {
            "gene": gene,
            "label": gene,
            "color": DARK,
            "marine_coef": np.nan,
            "aquatic_coef": np.nan,
            "display_coef": np.nan,
            "marine_marker": "",
            "aquatic_marker": "",
        }

    marine_coef = vals["marine_coef"]
    aquatic_coef = vals["aquatic_coef"]
    if context == "shared":
        display_coef = np.nanmean([marine_coef, aquatic_coef])
    elif context == "marine_specific":
        display_coef = marine_coef
    elif context == "aquatic_specific":
        display_coef = aquatic_coef
    else:
        display_coef = np.nan

    marine_marker = "#" if abs(marine_coef) > 0.1 else ""
    aquatic_marker = "*" if abs(aquatic_coef) > 0.1 else ""
    if context == "marine_specific":
        marker = marine_marker
    elif context == "aquatic_specific":
        marker = aquatic_marker
    else:
        marker = marine_marker + aquatic_marker

    return {
        "gene": gene,
        "label": gene + marker,
        "color": coef_color(display_coef) if np.isfinite(display_coef) else DARK,
        "marine_coef": marine_coef,
        "aquatic_coef": aquatic_coef,
        "display_coef": display_coef,
        "marine_marker": marine_marker,
        "aquatic_marker": aquatic_marker,
    }


def rich_gene_label(text: object, context: str, coef_lookup: dict[str, dict[str, float]]) -> list[dict[str, object]]:
    return [format_gene_label(gene, context, coef_lookup) for gene in split_genes(text)]


def label_text(parts: list[dict[str, object]]) -> str:
    return ", ".join(str(part["label"]) for part in parts)


def bubble_size(n: int) -> float:
    if n <= 0:
        return 0
    return 170 + (n ** 1.45) * 72


def build_plot_table() -> pd.DataFrame:
    module = pd.read_csv(MODULE_FILE, sep="\t")
    gene_level = pd.read_csv(GENE_LEVEL_FILE, sep="\t")

    module["module"] = module["module"].map(normalize_module_name)

    required_module_cols = [
        "module", "n_shared", "n_marine_specific", "n_aquatic_specific",
        "shared_genes", "marine_specific_genes", "aquatic_specific_genes",
    ]
    require_columns(module, required_module_cols, str(MODULE_FILE))
    require_columns(gene_level, ["gene", "lasso_group", "candidate_module", "marine_coef", "aquatic_coef"], str(GENE_LEVEL_FILE))
    gene_level["marine_coef"] = pd.to_numeric(gene_level["marine_coef"], errors="coerce").fillna(0.0)
    gene_level["aquatic_coef"] = pd.to_numeric(gene_level["aquatic_coef"], errors="coerce").fillna(0.0)
    coef_lookup = gene_level.set_index("gene")[["marine_coef", "aquatic_coef"]].to_dict(orient="index")

    missing_modules = [m for m in MODULE_ORDER if m not in set(module["module"])]
    if missing_modules:
        raise ValueError("Missing expected modules: " + "; ".join(missing_modules))

    module = module.set_index("module").loc[MODULE_ORDER].reset_index()

    rows = []
    for y_idx, row in module.iterrows():
        y = len(MODULE_ORDER) - y_idx
        for category, count_col, gene_col, color, x in CATEGORIES:
            n = int(row[count_col])
            if count_col == "n_shared":
                context = "shared"
            elif count_col == "n_marine_specific":
                context = "marine_specific"
            else:
                context = "aquatic_specific"
            gene_parts = rich_gene_label(row.get(gene_col, ""), context, coef_lookup)
            rows.append({
                "module": row["module"],
                "module_label": MODULE_LABELS[row["module"]],
                "gene_set_category": category,
                "gene_count": n,
                "gene_names": full_gene_label(row.get(gene_col, "")),
                "gene_label_with_markers": label_text(gene_parts),
                "gene_label_parts": gene_parts,
                "x_position": x,
                "y_position": y,
                "bubble_size": bubble_size(n),
                "color": color,
            })

    plot_table = pd.DataFrame(rows)
    plot_table.to_csv(PLOT_TABLE_FILE, sep="\t", index=False)

    print("Figure 4C module counts used")
    export_table = plot_table.copy()
    export_table["gene_label_parts"] = export_table["gene_label_parts"].map(
        lambda parts: ";".join(
            f"{part['gene']}|{part['label']}|{part['color']}|marine={part['marine_coef']}|aquatic={part['aquatic_coef']}|display={part['display_coef']}"
            for part in parts
        )
    )
    export_table.to_csv(PLOT_TABLE_FILE, sep="\t", index=False)

    print("Figure 4C module counts used")
    print(export_table[["module", "gene_set_category", "gene_count", "gene_label_with_markers"]].to_string(index=False))
    return plot_table


def add_colored_gene_label(ax, x: float, y: float, parts: list[dict[str, object]], fontsize: float = 6.2) -> None:
    if not parts:
        return
    children = []
    for i, part in enumerate(parts):
        text = str(part["label"]) + (", " if i < len(parts) - 1 else "")
        children.append(
            TextArea(
                text,
                textprops={
                    "fontsize": fontsize,
                    "color": part["color"],
                    "fontstyle": "italic",
                    "fontfamily": "DejaVu Sans",
                },
            )
        )
    packed = HPacker(children=children, align="center", pad=0, sep=0)
    anchored = AnchoredOffsetbox(
        loc="center left",
        child=packed,
        pad=0,
        frameon=False,
        bbox_to_anchor=(x, y),
        bbox_transform=ax.transData,
        borderpad=0,
    )
    anchored.set_zorder(5)
    anchored.set_clip_on(False)
    ax.add_artist(anchored)


def draw_panel(plot_table: pd.DataFrame, with_examples: bool, pdf_file: Path, svg_file: Path, png_file: Path | None = None) -> None:
    fig_width = 18.4 if with_examples else 14.8
    fig, ax = plt.subplots(figsize=(fig_width, 7.2))
    ax.set_facecolor("white")

    # Light row bands and vertical guides.
    for y in range(1, len(MODULE_ORDER) + 1):
        ax.axhline(y, color="#F0F0F0", lw=0.8, zorder=0)
    for x in [0.0, 1.28, 2.65]:
        ax.axvline(x, color=GRID, lw=0.9, zorder=0)

    # Module labels.
    for i, module_name in enumerate(MODULE_ORDER):
        y = len(MODULE_ORDER) - i
        ax.text(-0.86, y, MODULE_LABELS[module_name], ha="right", va="center",
                fontsize=9.3, color=DARK, linespacing=1.08)

    # Column headers.
    headers = [
        ("Shared aquatic foundation", SHARED_GREY),
        ("Marine-specific predictors", MARINE_BLUE),
        ("Aquatic-specific predictors", AQUATIC_TEAL),
    ]
    for x, (label, color) in zip([0.0, 1.28, 2.65], headers):
        header_label = label if not with_examples else label.replace(" predictors", "\npredictors").replace(" foundation", "\nfoundation")
        ax.text(x, len(MODULE_ORDER) + 0.38, header_label, ha="center", va="bottom",
                fontsize=9.2 if with_examples else 10.5, fontweight="bold",
                color=color, linespacing=0.9)

    # Bubbles and count labels.
    nonzero = plot_table[plot_table["gene_count"] > 0].copy()
    zero = plot_table[plot_table["gene_count"] == 0].copy()
    ax.scatter(
        nonzero["x_position"], nonzero["y_position"],
        s=nonzero["bubble_size"], c=nonzero["color"],
        edgecolors="#333333", linewidths=0.45, alpha=0.92, zorder=3,
    )
    for _, row in nonzero.iterrows():
        ax.text(row["x_position"], row["y_position"], str(row["gene_count"]),
                ha="center", va="center", fontsize=10.2, fontweight="bold", color=DARK, zorder=4)

    for _, row in zero.iterrows():
        ax.text(row["x_position"], row["y_position"], "0",
                ha="center", va="center", fontsize=8.0, color="#B0B0B0", zorder=4)

    if with_examples:
        for _, row in nonzero.iterrows():
            parts = row["gene_label_parts"]
            if not parts:
                continue
            # Place example labels to the right within each column block.
            add_colored_gene_label(ax, row["x_position"] + 0.18, row["y_position"], parts, fontsize=6.2)

    # Legend for bubble size.
    legend_counts = [1, 2, 4, 6, 9]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="#333333", markersize=np.sqrt(bubble_size(n)) / 3.4,
               label=str(n))
        for n in legend_counts
    ]
    ax.legend(handles=handles, title="# genes", frameon=False, loc="center left",
              bbox_to_anchor=(1.02, 0.44), fontsize=8.5, title_fontsize=9.5,
              labelspacing=1.15, handletextpad=1.0)

    title = "Module-annotated partition of selected predictors"
    subtitle = "Count-only matrix" if not with_examples else "All module-annotated gene names; red = positive/fast, blue = negative/slow; # marine |coef| > 0.1, * aquatic |coef| > 0.1"
    ax.text(-0.86, len(MODULE_ORDER) + 1.62, title, ha="left", va="bottom",
            fontsize=13.0, fontweight="bold", color=DARK)
    ax.text(-0.86, len(MODULE_ORDER) + 1.34, subtitle, ha="left", va="bottom",
            fontsize=9.3, color="#555555")

    ax.set_xlim(-1.05, 5.25 if with_examples else 3.05)
    ax.set_ylim(0.35, len(MODULE_ORDER) + 2.05)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.24 if with_examples else 0.30,
                        right=0.91 if with_examples else 0.88,
                        top=0.92, bottom=0.06)
    fig.savefig(pdf_file, transparent=True)
    fig.savefig(svg_file, transparent=True)
    if png_file is not None:
        fig.savefig(png_file, transparent=True, dpi=300)
    plt.close(fig)
    print(f"Saved: {pdf_file}")
    print(f"Saved: {svg_file}")
    if png_file is not None:
        print(f"Saved: {png_file}")


def main() -> None:
    plot_table = build_plot_table()
    draw_panel(plot_table, with_examples=False, pdf_file=COUNT_ONLY_PDF, svg_file=COUNT_ONLY_SVG, png_file=COUNT_ONLY_PNG)
    draw_panel(plot_table, with_examples=True, pdf_file=EXAMPLES_PDF, svg_file=EXAMPLES_SVG, png_file=EXAMPLES_PNG)
    draw_panel(plot_table, with_examples=True, pdf_file=FINAL_PDF, svg_file=FINAL_SVG, png_file=FINAL_PNG)
    print(f"Saved: {PLOT_TABLE_FILE}")


if __name__ == "__main__":
    main()
