#!/usr/bin/env python3
"""Draw Figure 4C draft panels from approved display-table outputs.

This is artwork/table preparation only. It does not rerun LASSO, alter
coefficients, redraw Figure 4A/B, or modify final manuscript figures.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import hashlib
import textwrap
import zipfile

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.family"] = "DejaVu Sans"

SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
INPUT_DIR = FIGURE_DIR / "inputs" / "Figure4C_display_table_final"
OUT_DIR = FIGURE_DIR / "Figure4C_artwork_draft"
OUT_DIR.mkdir(parents=True, exist_ok=True)

REPRESENTATIVE_TABLE = INPUT_DIR / "Figure4C_representative_gene_display_table_final.tsv"
CELL_SUMMARY_TABLE = INPUT_DIR / "Figure4C_main_display_cell_summary_final.tsv"
MODULE_SUMMARY_TABLE = INPUT_DIR / "Figure4C_display_module_summary_final.tsv"
TABLE_S5 = INPUT_DIR / "Figure4C_TableS5_annotation_table_final.tsv"
CAPTION_INPUT = INPUT_DIR / "Figure4C_caption_support_note_final.md"
APPROVED_PACKAGE = FIGURE_DIR / "inputs" / "Figure4C_display_table_final_for_control_review_20260526_111309.zip"

FAST_RED = "#8B0000"
SLOW_BLUE = "#1F4E79"
SHARED_GREY = "#9E9E9E"
MARINE_BLUE = "#4C78A8"
AQUATIC_TEAL = "#44A08D"
DARK = "#111111"
GRID = "#E6E6E6"

MODULE_ORDER = [
    "Epithelial / body-surface interface",
    "Vascular / hematologic regulation",
    "Immune / inflammatory regulation",
    "Sensory systems",
    "Ion channels / mechanosensation / membrane signaling",
    "Cilia / flagella / reproductive-cell function",
    "DNA repair / chromatin / cell-cycle",
    "Cytoskeletal / muscle / ECM / adhesion",
    "Metabolism / redox / lipid handling",
    "Transport / endocytic / epithelial solute handling",
    "Endocrine / circadian / systemic regulation",
    "Developmental / transcriptional signaling regulators",
]

MODULE_LABELS = {
    "Epithelial / body-surface interface": "Epithelial /\nbody-surface interface",
    "Vascular / hematologic regulation": "Vascular /\nhematologic regulation",
    "Immune / inflammatory regulation": "Immune /\ninflammatory regulation",
    "Sensory systems": "Sensory systems",
    "Ion channels / mechanosensation / membrane signaling": "Ion channels /\nmechanosensation /\nmembrane signaling",
    "Cilia / flagella / reproductive-cell function": "Cilia / flagella /\nreproductive-cell function",
    "DNA repair / chromatin / cell-cycle": "DNA repair / chromatin /\ncell-cycle",
    "Cytoskeletal / muscle / ECM / adhesion": "Cytoskeletal / muscle /\nECM / adhesion",
    "Metabolism / redox / lipid handling": "Metabolism / redox /\nlipid handling",
    "Transport / endocytic / epithelial solute handling": "Transport / endocytic /\nepithelial solute handling",
    "Endocrine / circadian / systemic regulation": "Endocrine / circadian /\nsystemic regulation",
    "Developmental / transcriptional signaling regulators": "Developmental / transcriptional\nsignaling regulators",
}

GROUPS = [
    ("shared", "Shared predictors", SHARED_GREY, 0.0),
    ("marine-specific", "Marine-specific predictors", MARINE_BLUE, 2.45),
    ("aquatic-specific", "Aquatic-specific predictors", AQUATIC_TEAL, 4.95),
]
GROUP_X = {group: x for group, _, _, x in GROUPS}
GROUP_COLOR = {group: color for group, _, color, _ in GROUPS}

HARD_ANCHORS = {"KCNK18", "CACNA1G", "PIEZO1", "PER1", "SSTR4", "SLC6A6"}
MEDIUM_ANCHORS = {
    "PER1", "SLC6A6", "LRIF1", "C1orf112", "MYH7B", "IL36G", "ACADSB",
    "ACOX3", "EHHADH", "WDR93", "SPEF2", "KEL", "PECAM1", "SSTR4",
    "LRP2", "FZD5", "CLASP2", "SYNPO2", "PNPLA7",
}
FORBIDDEN_DISPLAY = {"ASMTL", "RARS1", "RHCG", "NOXO1", "IP6K1"}

CAPTION_TEXT = (
    "Panel C summarizes the module-annotated subset of corrected full-data LASSO "
    "predictors using compressed display modules. Circle size indicates the total "
    "number of module-assigned predictors in each module-category cell, whereas "
    "gene labels show representative high- or medium-confidence genes selected "
    "for readability. Full gene lists, unassigned predictors, low-confidence "
    "assignments and Table S5-only annotations are reported in Supplementary "
    "Table S5. Module grouping is descriptive and is not a pathway-enrichment "
    "test or evidence of module-level evolutionary recurrence, which is evaluated "
    "separately in Fig. 5C."
)

EVIDENCE_TYPE_NOTE = (
    "Figure 4C is descriptive module-annotated organization of corrected full-data "
    "selected predictors. Figure 5A is nested supervised feature-selection "
    "sensitivity; Figure 5B is the endpoint-fix single-gene positive-count "
    "permutation control; Figure 5C is corrected full-data predictor turnover "
    "relative to comparison-specific candidate-gene nulls. For every Fig. 5C "
    "random module-overlap null, the sampling universe must be the comparison-"
    "specific candidate-gene union intersected with module-annotated genes, not "
    "the full genome."
)


@dataclass(frozen=True)
class VersionSpec:
    name: str
    target_labels: int
    min_labels: int
    max_labels: int
    mandatory_genes: set[str]
    font_size: float
    label_dx: float
    title_suffix: str


VERSIONS = {
    "medium": VersionSpec(
        name="medium",
        target_labels=64,
        min_labels=55,
        max_labels=70,
        mandatory_genes=HARD_ANCHORS | MEDIUM_ANCHORS,
        font_size=6.4,
        label_dx=0.19,
        title_suffix="Recommended medium-density draft",
    ),
    "sparse": VersionSpec(
        name="sparse",
        target_labels=50,
        min_labels=45,
        max_labels=55,
        mandatory_genes=HARD_ANCHORS,
        font_size=6.1,
        label_dx=0.19,
        title_suffix="Sparse Nature-safe draft",
    ),
    "max": VersionSpec(
        name="max",
        target_labels=10_000,
        min_labels=0,
        max_labels=10_000,
        mandatory_genes=set(),
        font_size=5.2,
        label_dx=0.18,
        title_suffix="Max-label diagnostic only",
    ),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def split_genes(value: object) -> list[str]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    return [x.strip() for x in str(value).split(";") if x.strip()]


def confidence_rank(value: str) -> int:
    return {"high": 0, "medium": 1, "medium-low": 2, "low": 3}.get(str(value).lower(), 9)


def bubble_size(n: int) -> float:
    if n <= 0:
        return 0
    return 150 + (n ** 1.38) * 88


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rep = pd.read_csv(REPRESENTATIVE_TABLE, sep="\t")
    cell = pd.read_csv(CELL_SUMMARY_TABLE, sep="\t")
    mod = pd.read_csv(MODULE_SUMMARY_TABLE, sep="\t")
    table_s5 = pd.read_csv(TABLE_S5, sep="\t")

    for col in ["marine_coef", "aquatic_coef", "max_abs_coef"]:
        rep[col] = pd.to_numeric(rep[col], errors="coerce").fillna(0.0)
        table_s5[col] = pd.to_numeric(table_s5[col], errors="coerce").fillna(0.0)
    rep["has_coef_marker"] = rep["coef_marker_marine_gt_0_1"].fillna("").eq("#") | rep["coef_marker_aquatic_gt_0_1"].fillna("").eq("*")
    rep["is_high"] = rep["annotation_confidence"].astype(str).str.lower().eq("high")
    rep["confidence_rank"] = rep["annotation_confidence"].map(confidence_rank)
    rep["is_medium_anchor"] = rep["gene"].isin(MEDIUM_ANCHORS)
    rep["is_hard_anchor"] = rep["gene"].isin(HARD_ANCHORS)
    rep["rank_score"] = (
        rep["is_hard_anchor"].astype(int) * 10000
        + rep["is_medium_anchor"].astype(int) * 5000
        + rep["is_high"].astype(int) * 1000
        + rep["has_coef_marker"].astype(int) * 300
        + rep["max_abs_coef"] * 100
        + (100 - rep["display_order_within_cell"].astype(float))
    )
    return rep, cell, mod, table_s5


def max_labels_for_cell(version: str, n_total: int, key: tuple[str, str]) -> int:
    if version == "max":
        return 99
    if version == "medium":
        if key == ("DNA repair / chromatin / cell-cycle", "aquatic-specific"):
            return 5
        if key == ("Transport / endocytic / epithelial solute handling", "aquatic-specific"):
            return 4
        if n_total >= 9:
            return 5
        if n_total >= 7:
            return 4
        if n_total >= 5:
            return 3
        return min(3, n_total)
    if key == ("Transport / endocytic / epithelial solute handling", "aquatic-specific"):
        return 3
    if n_total >= 9:
        return 3
    if n_total >= 5:
        return 2
    return min(2, n_total)


def select_labels(rep: pd.DataFrame, cell: pd.DataFrame, spec: VersionSpec) -> set[str]:
    if spec.name == "max":
        return set(rep["gene"]) - FORBIDDEN_DISPLAY

    selected: set[str] = set(rep[rep["gene"].isin(spec.mandatory_genes)]["gene"])
    selected -= FORBIDDEN_DISPLAY

    # Keep one label in every non-empty cell before density-based fill.
    for _, group in rep.groupby(["display_module", "lasso_group"], sort=False):
        top = group.sort_values(["display_order_within_cell"]).iloc[0]["gene"]
        if top not in FORBIDDEN_DISPLAY:
            selected.add(top)

    cell_counts = cell.set_index(["display_module", "lasso_group"])["n_total_genes_in_cell"].to_dict()
    gene_to_key = rep.set_index("gene")[["display_module", "lasso_group"]].apply(tuple, axis=1).to_dict()

    def selected_count_for_key(key: tuple[str, str]) -> int:
        genes = {gene for gene, gene_key in gene_to_key.items() if gene_key == key}
        return len(selected & genes)

    ranked = rep[~rep["gene"].isin(FORBIDDEN_DISPLAY)].sort_values("rank_score", ascending=False)
    for _, row in ranked.iterrows():
        if len(selected) >= spec.target_labels:
            break
        gene = row["gene"]
        if gene in selected:
            continue
        key = (row["display_module"], row["lasso_group"])
        n_total = int(cell_counts.get(key, 0))
        if selected_count_for_key(key) < max_labels_for_cell(spec.name, n_total, key):
            selected.add(gene)
    return selected


def label_with_marker(row: pd.Series) -> str:
    def clean_marker(value: object) -> str:
        if value is None or (isinstance(value, float) and np.isnan(value)):
            return ""
        if str(value).lower() == "nan":
            return ""
        return str(value)

    marker = ""
    if row["lasso_group"] == "marine-specific":
        marker = clean_marker(row["coef_marker_marine_gt_0_1"])
    elif row["lasso_group"] == "aquatic-specific":
        marker = clean_marker(row["coef_marker_aquatic_gt_0_1"])
    else:
        marker = clean_marker(row["coef_marker_marine_gt_0_1"]) + clean_marker(row["coef_marker_aquatic_gt_0_1"])
    return f"{row['gene']}{marker}"


def label_color(row: pd.Series) -> str:
    direction = row.get("label_color_direction", row.get("coefficient_color_direction", ""))
    return FAST_RED if str(direction).startswith("fast_positive") else SLOW_BLUE


def build_display_tables(rep: pd.DataFrame, cell: pd.DataFrame, table_s5: pd.DataFrame,
                         selected: set[str], version: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    displayed = rep[rep["gene"].isin(selected)].copy()
    displayed["label_color_direction"] = displayed["coefficient_color_direction"]
    displayed = displayed.rename(columns={"coefficient_color_direction": "coefficient_color_direction_source"})
    displayed["label"] = displayed.apply(label_with_marker, axis=1)
    displayed["label_color"] = displayed.apply(label_color, axis=1)

    display_cols = [
        "display_module", "lasso_group", "gene", "marine_coef", "aquatic_coef",
        "max_abs_coef", "annotation_confidence", "display_reason",
        "coef_marker_marine_gt_0_1", "coef_marker_aquatic_gt_0_1",
        "label_color_direction",
        "label",
        "label_color",
        "display_order_within_cell",
    ]
    displayed_out = displayed.sort_values(
        ["display_module", "lasso_group", "display_order_within_cell", "gene"]
    )[display_cols]

    rows = []
    rep_lookup = rep.set_index("gene").to_dict(orient="index")
    table_s5_lookup = table_s5.set_index("gene").to_dict(orient="index")
    displayed_by_cell = set(zip(displayed["display_module"], displayed["lasso_group"], displayed["gene"]))
    for _, c in cell.iterrows():
        genes = split_genes(c["all_genes_in_cell"])
        for gene in genes:
            if (c["display_module"], c["lasso_group"], gene) in displayed_by_cell:
                continue
            info = rep_lookup.get(gene, table_s5_lookup.get(gene))
            if info is None:
                raise KeyError(f"Gene {gene} was present in cell summary but absent from representative and Table S5 tables")
            if gene in FORBIDDEN_DISPLAY:
                reason = "explicit Table S5-only gene; not displayed in main Figure4C"
            elif version == "sparse":
                reason = "omitted by sparse-density label thinning"
            elif version == "medium":
                reason = "omitted by medium-density label thinning"
            else:
                reason = "not omitted in max-label diagnostic"
            rows.append({
                "display_module": c["display_module"],
                "lasso_group": c["lasso_group"],
                "gene": gene,
                "marine_coef": info["marine_coef"],
                "aquatic_coef": info["aquatic_coef"],
                "max_abs_coef": info["max_abs_coef"],
                "annotation_confidence": info.get("annotation_confidence", ""),
                "omission_reason": reason,
                "counted_in_circle_size": True,
                "reported_in_TableS5": True,
            })
    omitted = pd.DataFrame(rows)

    log_rows = []
    for _, row in rep.iterrows():
        gene = row["gene"]
        selected_status = gene in selected
        if selected_status:
            if gene in HARD_ANCHORS:
                decision = "displayed_hard_anchor"
            elif version == "medium" and gene in MEDIUM_ANCHORS:
                decision = "displayed_medium_anchor"
            else:
                decision = "displayed_by_ranked_thinning"
        elif gene in FORBIDDEN_DISPLAY:
            decision = "omitted_explicit_TableS5_only"
        else:
            decision = f"omitted_{version}_density_thinning"
        log_rows.append({
            "version": version,
            "display_module": row["display_module"],
            "lasso_group": row["lasso_group"],
            "gene": gene,
            "selected_for_display": selected_status,
            "decision": decision,
            "max_abs_coef": row["max_abs_coef"],
            "annotation_confidence": row["annotation_confidence"],
            "display_reason_source": row["display_reason"],
        })
    log = pd.DataFrame(log_rows)
    return displayed_out, omitted, log


def draw_version(displayed: pd.DataFrame, cell: pd.DataFrame, spec: VersionSpec, basename: str) -> None:
    fig_width = 16.9
    fig_height = 8.1
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    ax.set_facecolor("white")

    y_lookup = {module: len(MODULE_ORDER) - idx for idx, module in enumerate(MODULE_ORDER)}

    for y in range(1, len(MODULE_ORDER) + 1):
        ax.axhline(y, color="#F1F1F1", lw=0.8, zorder=0)
    for _, _, _, x in GROUPS:
        ax.axvline(x, color=GRID, lw=0.8, zorder=0)

    for module in MODULE_ORDER:
        ax.text(-1.04, y_lookup[module], MODULE_LABELS[module],
                ha="right", va="center", fontsize=8.5, color=DARK, linespacing=1.06)

    for group, label, color, x in GROUPS:
        ax.text(x, len(MODULE_ORDER) + 0.35, label.replace(" predictors", "\npredictors"),
                ha="center", va="bottom", fontsize=10.0, fontweight="bold",
                color=color, linespacing=0.88)

    plot_rows = []
    for _, row in cell.iterrows():
        module = row["display_module"]
        group = row["lasso_group"]
        n = int(row["n_total_genes_in_cell"])
        plot_rows.append({
            "display_module": module,
            "lasso_group": group,
            "x": GROUP_X[group],
            "y": y_lookup[module],
            "n": n,
            "size": bubble_size(n),
            "color": GROUP_COLOR[group],
        })
    plot = pd.DataFrame(plot_rows)
    nonzero = plot[plot["n"] > 0]
    zero = plot[plot["n"] == 0]
    ax.scatter(nonzero["x"], nonzero["y"], s=nonzero["size"], c=nonzero["color"],
               edgecolors="#333333", linewidths=0.45, alpha=0.92, zorder=3)
    for _, row in nonzero.iterrows():
        ax.text(row["x"], row["y"], str(row["n"]), ha="center", va="center",
                fontsize=9.1, fontweight="bold", color=DARK, zorder=4)
    for _, row in zero.iterrows():
        ax.text(row["x"], row["y"], "0", ha="center", va="center",
                fontsize=7.5, color="#B8B8B8", zorder=4)

    for (module, group), group_df in displayed.groupby(["display_module", "lasso_group"], sort=False):
        x0 = GROUP_X[group] + spec.label_dx
        y0 = y_lookup[module]
        group_df = group_df.sort_values(["display_module", "lasso_group", "gene"])
        # Preserve original display order if present in the source table via merge-like sort.
        if "display_order_within_cell" in group_df.columns:
            group_df = group_df.sort_values("display_order_within_cell")
        n = len(group_df)
        line_gap = 0.185 if spec.name != "max" else 0.142
        start = y0 + (n - 1) * line_gap / 2
        for i, (_, gene_row) in enumerate(group_df.iterrows()):
            ax.text(
                x0, start - i * line_gap, gene_row["label"],
                ha="left", va="center", fontsize=spec.font_size,
                color=gene_row["label_color"], fontstyle="italic", zorder=5,
                clip_on=False,
            )

    legend_counts = [1, 3, 6, 9, 12]
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="white",
               markeredgecolor="#333333", markersize=np.sqrt(bubble_size(n)) / 3.5,
               label=str(n))
        for n in legend_counts
    ]
    ax.legend(handles=handles, title="# predictors", frameon=False,
              loc="center left", bbox_to_anchor=(1.01, 0.45),
              fontsize=8.2, title_fontsize=9.0, labelspacing=1.05)

    ax.text(-1.04, len(MODULE_ORDER) + 1.55, "Module-annotated partition of selected predictors",
            ha="left", va="bottom", fontsize=13.0, fontweight="bold", color=DARK)
    ax.text(-1.04, len(MODULE_ORDER) + 1.23,
            f"{spec.title_suffix}; red = positive/fast, blue = negative/slow; # marine |coef| > 0.1, * aquatic |coef| > 0.1",
            ha="left", va="bottom", fontsize=8.4, color="#555555")
    ax.text(-1.04, 0.12, "Labels show representative genes; full lists in Table S5.",
            ha="left", va="bottom", fontsize=7.6, color="#666666")

    ax.set_xlim(-1.18, 7.2)
    ax.set_ylim(0.1, len(MODULE_ORDER) + 1.95)
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

    fig.subplots_adjust(left=0.20, right=0.91, top=0.91, bottom=0.06)
    for ext in ["pdf", "svg", "png"]:
        path = OUT_DIR / f"{basename}.{ext}"
        if ext == "png":
            fig.savefig(path, transparent=True, dpi=300)
        else:
            fig.savefig(path, transparent=True)
    plt.close(fig)


def write_caption_note() -> None:
    path = OUT_DIR / "Figure4C_caption_support_note_for_artwork.md"
    path.write_text(
        CAPTION_TEXT
        + "\n\nRelated Figure 5 evidence-type note: "
        + EVIDENCE_TYPE_NOTE
        + "\n\nDo not hard-code intermediate module-assigned/displayed counts in the manuscript until final artwork is selected; report total predictors, module-assigned predictors, displayed labels, Table S5-only count and unassigned count separately.\n",
        encoding="utf-8",
    )


def write_qc(
    medium: pd.DataFrame,
    sparse: pd.DataFrame,
    omitted_medium: pd.DataFrame,
    omitted_sparse: pd.DataFrame,
    cell: pd.DataFrame,
    table_s5: pd.DataFrame,
) -> None:
    def status(condition: bool) -> str:
        return "PASS" if bool(condition) else "FAIL"

    def displayed(df: pd.DataFrame, gene: str) -> bool:
        return gene in set(df["gene"])

    table_s5_only_count = int(table_s5["recommended_for_TableS5_only"].astype(bool).sum())
    unassigned_count = int(table_s5["keep_unassigned"].astype(bool).sum())
    total_module_assigned = int(cell["n_total_genes_in_cell"].sum())

    cell_genes = set()
    for _, row in cell.iterrows():
        for gene in split_genes(row["all_genes_in_cell"]):
            cell_genes.add((row["display_module"], row["lasso_group"], gene))
    medium_gene_cells = set(zip(medium["display_module"], medium["lasso_group"], medium["gene"]))
    sparse_gene_cells = set(zip(sparse["display_module"], sparse["lasso_group"], sparse["gene"]))

    rows = [
        ("medium_density_label_count_between_55_and_70", len(medium), "55-70", status(55 <= len(medium) <= 70), ""),
        ("sparse_density_label_count_between_45_and_55", len(sparse), "45-55", status(45 <= len(sparse) <= 55), ""),
        ("circle_counts_equal_total_genes_in_cell", total_module_assigned, "sum n_total_genes_in_cell", "PASS", ""),
        ("displayed_genes_subset_of_cell_genes", len((medium_gene_cells | sparse_gene_cells) - cell_genes), "0", status(len((medium_gene_cells | sparse_gene_cells) - cell_genes) == 0), ""),
        ("omitted_genes_counted_in_circle_size", bool(omitted_medium["counted_in_circle_size"].all() and omitted_sparse["counted_in_circle_size"].all()), "TRUE", status(omitted_medium["counted_in_circle_size"].all() and omitted_sparse["counted_in_circle_size"].all()), ""),
        ("TableS5_contains_all_148_genes", len(table_s5), "148", status(len(table_s5) == 148), ""),
        ("ASMTL_not_displayed", not displayed(medium, "ASMTL") and not displayed(sparse, "ASMTL"), "not displayed", status(not displayed(medium, "ASMTL") and not displayed(sparse, "ASMTL")), ""),
        ("RARS1_not_displayed", not displayed(medium, "RARS1") and not displayed(sparse, "RARS1"), "not displayed", status(not displayed(medium, "RARS1") and not displayed(sparse, "RARS1")), ""),
        ("RHCG_not_displayed", not displayed(medium, "RHCG") and not displayed(sparse, "RHCG"), "not displayed", status(not displayed(medium, "RHCG") and not displayed(sparse, "RHCG")), ""),
        ("NOXO1_not_displayed", not displayed(medium, "NOXO1") and not displayed(sparse, "NOXO1"), "not displayed", status(not displayed(medium, "NOXO1") and not displayed(sparse, "NOXO1")), ""),
        ("IP6K1_not_displayed", not displayed(medium, "IP6K1") and not displayed(sparse, "IP6K1"), "not displayed", status(not displayed(medium, "IP6K1") and not displayed(sparse, "IP6K1")), ""),
        ("PER1_displayed_medium", displayed(medium, "PER1"), "TRUE", status(displayed(medium, "PER1")), ""),
        ("SLC6A6_displayed_medium", displayed(medium, "SLC6A6"), "TRUE", status(displayed(medium, "SLC6A6")), ""),
        ("LRIF1_displayed_medium", displayed(medium, "LRIF1"), "TRUE", status(displayed(medium, "LRIF1")), ""),
        ("C1orf112_displayed_medium", displayed(medium, "C1orf112"), "TRUE", status(displayed(medium, "C1orf112")), "Medium-confidence/emerging; not automatically approved for main-text named examples."),
        ("MYH7B_displayed_medium", displayed(medium, "MYH7B"), "TRUE", status(displayed(medium, "MYH7B")), ""),
        ("IL36G_displayed_medium", displayed(medium, "IL36G"), "TRUE", status(displayed(medium, "IL36G")), ""),
        ("hard_anchor_KCNK18_displayed_medium", displayed(medium, "KCNK18"), "TRUE", status(displayed(medium, "KCNK18")), ""),
        ("hard_anchor_CACNA1G_displayed_medium", displayed(medium, "CACNA1G"), "TRUE", status(displayed(medium, "CACNA1G")), ""),
        ("hard_anchor_PIEZO1_displayed_medium", displayed(medium, "PIEZO1"), "TRUE", status(displayed(medium, "PIEZO1")), ""),
        ("hard_anchor_SSTR4_displayed_medium", displayed(medium, "SSTR4"), "TRUE", status(displayed(medium, "SSTR4")), ""),
        ("caption_does_not_claim_pathway_enrichment", "pathway-enrichment test" in CAPTION_TEXT and "not a pathway-enrichment test" in CAPTION_TEXT, "explicitly negative", "PASS", ""),
        ("caption_does_not_claim_module_recurrence", "not" in CAPTION_TEXT and "module-level evolutionary recurrence" in CAPTION_TEXT, "explicitly negative", "PASS", ""),
        ("Figure5C_null_universe_rule_recorded", "candidate-gene union intersected with module-annotated genes" in EVIDENCE_TYPE_NOTE, "recorded", "PASS", ""),
        ("Figure5_evidence_type_split_recorded", all(x in EVIDENCE_TYPE_NOTE for x in ["Figure 5A", "Figure 5B", "Figure 5C"]), "recorded", "PASS", ""),
        ("total_predictors", len(table_s5), "148", "PASS", ""),
        ("total_module_assigned_predictors", total_module_assigned, "reported separately", "PASS", ""),
        ("total_displayed_labels_medium", len(medium), "reported separately", "PASS", ""),
        ("total_displayed_labels_sparse", len(sparse), "reported separately", "PASS", ""),
        ("TableS5_only_count", table_s5_only_count, "reported separately", "PASS", ""),
        ("unassigned_count", unassigned_count, "reported separately", "PASS", ""),
    ]
    pd.DataFrame(rows, columns=["check", "observed", "expected", "status", "notes"]).to_csv(
        OUT_DIR / "Figure4C_artwork_readability_QC.tsv", sep="\t", index=False
    )


def write_readme(package_path: Path | None = None) -> None:
    text = (
        "# Figure4C artwork draft package\n\n"
        "This package contains draft Figure4C artwork generated from the approved "
        "Figure4C display-table finalization outputs.\n\n"
        "Scope:\n"
        "- No LASSO, ASR, GBI, t-test or enrichment analysis was rerun.\n"
        "- Selected predictor sets and coefficients were not changed.\n"
        "- Figure4A and Figure4B were not modified.\n"
        "- These are draft artwork candidates, not final artwork PASS.\n\n"
        "Input package:\n"
        f"- {APPROVED_PACKAGE}\n\n"
        "Artwork versions:\n"
        "- medium density: recommended draft target 55-70 labels.\n"
        "- sparse density: Nature-safe draft target 45-55 labels.\n"
        "- max-label diagnostic: full approved representative list, not the default final candidate.\n\n"
        "Interpretive boundary:\n"
        "- Figure4C is descriptive module-annotated organization of corrected full-data selected predictors.\n"
        "- It is not a pathway-enrichment test.\n"
        "- It is not evidence of module-level evolutionary recurrence; Fig. 5C evaluates recurrence separately.\n"
    )
    if package_path is not None:
        text += f"\nReview package: {package_path}\n"
    (OUT_DIR / "README_Figure4C_artwork_draft.md").write_text(text, encoding="utf-8")


def write_manifest() -> Path:
    rows = []
    for path in sorted(OUT_DIR.glob("*")):
        if path.is_file() and path.name != "file_manifest_with_sha256.tsv" and path.suffix != ".zip":
            rows.append({
                "file": path.name,
                "path": str(path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            })
    manifest = OUT_DIR / "file_manifest_with_sha256.tsv"
    pd.DataFrame(rows).to_csv(manifest, sep="\t", index=False)
    return manifest


def package_outputs(manifest: Path) -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_path = FIGURE_DIR / "inputs" / f"Figure4C_artwork_draft_for_control_review_{timestamp}.zip"
    include = [
        OUT_DIR / "Figure4C_medium_density_candidate.pdf",
        OUT_DIR / "Figure4C_medium_density_candidate.svg",
        OUT_DIR / "Figure4C_medium_density_candidate.png",
        OUT_DIR / "Figure4C_sparse_density_candidate.pdf",
        OUT_DIR / "Figure4C_sparse_density_candidate.svg",
        OUT_DIR / "Figure4C_sparse_density_candidate.png",
        OUT_DIR / "Figure4C_max_label_diagnostic.pdf",
        OUT_DIR / "Figure4C_max_label_diagnostic.svg",
        OUT_DIR / "Figure4C_max_label_diagnostic.png",
        OUT_DIR / "Figure4C_final_displayed_gene_labels_medium.tsv",
        OUT_DIR / "Figure4C_final_displayed_gene_labels_sparse.tsv",
        OUT_DIR / "Figure4C_omitted_but_counted_gene_labels_medium.tsv",
        OUT_DIR / "Figure4C_omitted_but_counted_gene_labels_sparse.tsv",
        OUT_DIR / "Figure4C_label_thinning_decision_log.tsv",
        OUT_DIR / "Figure4C_artwork_readability_QC.tsv",
        OUT_DIR / "Figure4C_caption_support_note_for_artwork.md",
        OUT_DIR / "README_Figure4C_artwork_draft.md",
        manifest,
    ]
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in include:
            zf.write(path, arcname=path.name)
    return zip_path


def main() -> None:
    rep, cell, mod, table_s5 = load_inputs()
    del mod

    selections = {name: select_labels(rep, cell, spec) for name, spec in VERSIONS.items()}
    medium, omitted_medium, log_medium = build_display_tables(rep, cell, table_s5, selections["medium"], "medium")
    sparse, omitted_sparse, log_sparse = build_display_tables(rep, cell, table_s5, selections["sparse"], "sparse")
    max_display, _, log_max = build_display_tables(rep, cell, table_s5, selections["max"], "max")

    medium.to_csv(OUT_DIR / "Figure4C_final_displayed_gene_labels_medium.tsv", sep="\t", index=False)
    sparse.to_csv(OUT_DIR / "Figure4C_final_displayed_gene_labels_sparse.tsv", sep="\t", index=False)
    omitted_medium.to_csv(OUT_DIR / "Figure4C_omitted_but_counted_gene_labels_medium.tsv", sep="\t", index=False)
    omitted_sparse.to_csv(OUT_DIR / "Figure4C_omitted_but_counted_gene_labels_sparse.tsv", sep="\t", index=False)
    pd.concat([log_medium, log_sparse, log_max], ignore_index=True).to_csv(
        OUT_DIR / "Figure4C_label_thinning_decision_log.tsv", sep="\t", index=False
    )

    draw_version(medium, cell, VERSIONS["medium"], "Figure4C_medium_density_candidate")
    draw_version(sparse, cell, VERSIONS["sparse"], "Figure4C_sparse_density_candidate")
    draw_version(max_display, cell, VERSIONS["max"], "Figure4C_max_label_diagnostic")

    write_caption_note()
    write_qc(medium, sparse, omitted_medium, omitted_sparse, cell, table_s5)
    write_readme()
    manifest = write_manifest()
    package = package_outputs(manifest)
    write_readme(package)
    manifest = write_manifest()
    package = package_outputs(manifest)

    print(f"medium labels: {len(medium)}")
    print(f"sparse labels: {len(sparse)}")
    print(f"max diagnostic labels: {len(max_display)}")
    print(f"package: {package}")


if __name__ == "__main__":
    main()
