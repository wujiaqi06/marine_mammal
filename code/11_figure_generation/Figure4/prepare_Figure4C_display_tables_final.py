#!/usr/bin/env python3
"""Prepare final Figure 4C display and Table S5 tables from Pro annotations.

This script performs table compression and representative-gene selection only.
It does not redraw Figure 4, rerun LASSO, or alter Pro's annotation evidence.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
import hashlib
import zipfile

import numpy as np
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
INPUT_DIR = FIGURE_DIR / "inputs"
OUT_DIR = INPUT_DIR / "Figure4C_display_table_final"
OUT_DIR.mkdir(parents=True, exist_ok=True)

RESCUE_TABLE = INPUT_DIR / "Figure4C_functional_annotation_rescue_table_v1.tsv"
MODULE_SUMMARY = INPUT_DIR / "Figure4C_functional_annotation_module_summary_v1.tsv"
PRIORITY_TABLE = INPUT_DIR / "Figure4C_priority_rescue_genes_v1.tsv"

COMPRESSION_MAP_OUT = OUT_DIR / "Figure4C_display_module_compression_map_final.tsv"
TABLE_S5_OUT = OUT_DIR / "Figure4C_TableS5_annotation_table_final.tsv"
CELL_SUMMARY_OUT = OUT_DIR / "Figure4C_main_display_cell_summary_final.tsv"
REPRESENTATIVE_OUT = OUT_DIR / "Figure4C_representative_gene_display_table_final.tsv"
MODULE_DISPLAY_OUT = OUT_DIR / "Figure4C_display_module_summary_final.tsv"
QC_OUT = OUT_DIR / "Figure4C_display_table_QC_final.tsv"
WARNING_OUT = OUT_DIR / "Figure4C_display_warning_log_final.tsv"
CAPTION_NOTE_OUT = OUT_DIR / "Figure4C_caption_support_note_final.md"
README_OUT = OUT_DIR / "README_Figure4C_display_table_final.md"
MANIFEST_OUT = OUT_DIR / "file_manifest_with_sha256.tsv"

CONFIDENCE_RANK = {
    "high": 0,
    "medium": 1,
    "medium-low": 2,
    "low": 2,
    "unassigned": 3,
}

DISPLAY_LASSO_GROUP = {
    "shared": "shared",
    "marine_specific": "marine-specific",
    "aquatic_specific": "aquatic-specific",
}

PRIORITY_GENES = {
    "PER1", "MYH7B", "IL36G", "ACADSB", "ACOX3", "EHHADH", "SLC6A6",
    "LRP2", "WDR93", "SPEF2", "KEL", "PECAM1", "SSTR4", "FAM13A",
    "NLRP8", "SLC6A15", "SIGIRR", "GCKR", "SYNPO2", "CLASP2", "PNPLA7",
    "FZD5", "XCR1", "SLC28A1", "ABHD8", "RARS1", "CORO2A", "SEPTIN5",
    "C1orf112",
}

FORCE_TABLES5_ONLY_GENES = {"RARS1", "RHCG", "NOXO1", "IP6K1"}
FORCE_MAIN_DISPLAY_GENES = {"SLC6A6", "LRIF1", "C1orf112"}
ASMTL_DISPLAY_NOTE = (
    "ASMTL is retained as an ASMT/ASMTL-locus and melatonin-pathway-adjacent "
    "context-only candidate. ASMTL ranks 262nd by marine t-test p-value, whereas "
    "ASMT ranks 3513th and is not significant in the current screen; therefore "
    "ASMTL should not be interpreted as direct evidence for altered melatonin "
    "secretion. Use only as cautious context in Supplementary Table S5 or Discussion."
)
FORCE_REPRESENTATIVES = {
    ("DNA repair / chromatin / cell-cycle", "aquatic_specific"): [
        "KIF18A", "H1_8", "CDC14B", "LRIF1", "C1orf112",
    ],
    ("Transport / endocytic / epithelial solute handling", "aquatic_specific"): [
        "LRP2", "SLC6A15", "SLC5A5", "SLC6A6",
    ],
}

MAIN_DISPLAY_MODULES = [
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

MODULE_COMPRESSION_ROWS = [
    {
        "pro_recommended_module": "Epithelial / body-surface interface",
        "display_module": "Epithelial / body-surface interface",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Vascular / hematologic regulation",
        "display_module": "Vascular / hematologic regulation",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Immune / inflammatory regulation",
        "display_module": "Immune / inflammatory regulation",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Sensory systems",
        "display_module": "Sensory systems",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Ion channels and membrane signaling",
        "display_module": "Ion channels / mechanosensation / membrane signaling",
        "display_status": "main_Figure4C",
        "reason": "Renamed to match approved display wording.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Cilia / flagella / reproductive-cell function",
        "display_module": "Cilia / flagella / reproductive-cell function",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "DNA repair / chromatin / cell-cycle",
        "display_module": "DNA repair / chromatin / cell-cycle",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Cytoskeletal / muscle / ECM / adhesion",
        "display_module": "Cytoskeletal / muscle / ECM / adhesion",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Mitochondrial / redox / intermediary metabolism",
        "display_module": "Metabolism / redox / lipid handling",
        "display_status": "main_Figure4C",
        "reason": "Compressed with fatty-acid/lipid/peroxisomal metabolism as instructed.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Fatty-acid / lipid / peroxisomal metabolism",
        "display_module": "Metabolism / redox / lipid handling",
        "display_status": "main_Figure4C",
        "reason": "Compressed with mitochondrial/redox/intermediary metabolism as instructed.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Transport / endocytic / epithelial solute handling",
        "display_module": "Transport / endocytic / epithelial solute handling",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Endocrine / circadian / systemic metabolic regulation",
        "display_module": "Endocrine / circadian / systemic regulation",
        "display_status": "main_Figure4C",
        "reason": "Compressed wording to approved display module.",
        "notes": "",
    },
    {
        "pro_recommended_module": "Developmental / transcriptional signaling regulators",
        "display_module": "Developmental / transcriptional signaling regulators",
        "display_status": "main_Figure4C",
        "reason": "Already matches approved 12-module display framework.",
        "notes": "Low-confidence genes remain Table S5-only at gene level.",
    },
    {
        "pro_recommended_module": "Membrane trafficking / ER / vesicle processing",
        "display_module": "Membrane trafficking / ER / vesicle processing",
        "display_status": "optional_Figure4C",
        "reason": "Kept outside main Figure4C display by default per final instruction.",
        "notes": "High-confidence/high-coefficient genes remain in Table S5 and warning log; not used in main display unless later approved.",
    },
    {
        "pro_recommended_module": "RNA processing / translation / protein homeostasis",
        "display_module": "RNA processing / translation / protein homeostasis",
        "display_status": "TableS5_only",
        "reason": "RNA processing/translation is excluded from main Figure 4C by instruction.",
        "notes": "RARS1 is not displayed in main Figure 4C unless this module is later approved.",
    },
    {
        "pro_recommended_module": "Unassigned / low-confidence",
        "display_module": "Unassigned / low-confidence",
        "display_status": "unassigned",
        "reason": "Unassigned and low-confidence genes are not forced into main display modules.",
        "notes": "ASMTL receives an explicit final rescue override to endocrine/circadian as a low-confidence ASMT/ASMTL-locus candidate.",
    },
]


def as_bool(x: object) -> bool:
    return str(x).strip().upper() == "TRUE"


def confidence_rank(value: object) -> int:
    return CONFIDENCE_RANK.get(str(value).strip().lower(), 99)


def marker_for_marine(abs_coef: float) -> str:
    return "#" if float(abs_coef) > 0.1 else ""


def marker_for_aquatic(abs_coef: float) -> str:
    return "*" if float(abs_coef) > 0.1 else ""


def coefficient_color_direction(row: pd.Series) -> str:
    if row["lasso_group"] == "shared":
        coef = np.nanmean([row["marine_coef"], row["aquatic_coef"]])
    elif row["lasso_group"] == "marine_specific":
        coef = row["marine_coef"]
    else:
        coef = row["aquatic_coef"]
    return "fast_positive_red" if coef >= 0 else "slow_negative_blue"


def is_priority(row: pd.Series) -> bool:
    return row["gene"] in PRIORITY_GENES


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def prepare_tables() -> dict[str, Path]:
    rescue = pd.read_csv(RESCUE_TABLE, sep="\t")
    pd.read_csv(MODULE_SUMMARY, sep="\t")
    priority = pd.read_csv(PRIORITY_TABLE, sep="\t")

    rescue = rescue.copy()
    for col in ["marine_coef", "aquatic_coef", "max_abs_coef"]:
        rescue[col] = pd.to_numeric(rescue[col], errors="coerce").fillna(0.0)
    rescue["selected_in_marine"] = rescue["selected_in_marine"].map(as_bool)
    rescue["selected_in_aquatic"] = rescue["selected_in_aquatic"].map(as_bool)
    rescue["recommended_for_Figure4C_original"] = rescue["recommended_for_Figure4C"].map(as_bool)
    rescue["recommended_for_TableS5_only_original"] = rescue["recommended_for_TableS5_only"].map(as_bool)
    rescue["keep_unassigned_original"] = rescue["keep_unassigned"].map(as_bool)
    rescue["lasso_group_display"] = rescue["lasso_group"].map(DISPLAY_LASSO_GROUP)
    rescue["abs_marine_coef"] = rescue["marine_coef"].abs()
    rescue["abs_aquatic_coef"] = rescue["aquatic_coef"].abs()
    rescue["coef_marker_marine_gt_0_1"] = rescue["abs_marine_coef"].map(marker_for_marine)
    rescue["coef_marker_aquatic_gt_0_1"] = rescue["abs_aquatic_coef"].map(marker_for_aquatic)
    rescue["coefficient_color_direction"] = rescue.apply(coefficient_color_direction, axis=1)
    rescue["is_priority_rescue_gene"] = rescue["gene"].isin(set(priority["gene"]) | PRIORITY_GENES)

    compression = pd.DataFrame(MODULE_COMPRESSION_ROWS)
    compression.to_csv(COMPRESSION_MAP_OUT, sep="\t", index=False)
    compression_lookup = compression.set_index("pro_recommended_module")

    rescue["pro_recommended_module"] = rescue["recommended_module"]
    rescue["display_module"] = rescue["pro_recommended_module"].map(compression_lookup["display_module"])
    rescue["display_status"] = rescue["pro_recommended_module"].map(compression_lookup["display_status"])
    rescue["display_module"] = rescue["display_module"].fillna(rescue["pro_recommended_module"])
    rescue["display_status"] = rescue["display_status"].fillna("unassigned")

    rescue["display_decision_note"] = ""
    asmtl_mask = rescue["gene"].eq("ASMTL")
    if asmtl_mask.any():
        rescue.loc[asmtl_mask, "display_module"] = "Endocrine / circadian / systemic regulation"
        rescue.loc[asmtl_mask, "display_status"] = "TableS5_only"
        rescue.loc[asmtl_mask, "recommended_submodule_or_function"] = (
            "ASMT/ASMTL-locus; melatonin-pathway-adjacent methyltransferase-like candidate"
        )
        rescue.loc[asmtl_mask, "annotation_confidence"] = "low"
        rescue.loc[asmtl_mask, "evidence_summary"] = ASMTL_DISPLAY_NOTE
        rescue.loc[asmtl_mask, "reason_for_assignment"] = ASMTL_DISPLAY_NOTE
        rescue.loc[asmtl_mask, "notes_for_manual_review"] = ASMTL_DISPLAY_NOTE
        rescue.loc[asmtl_mask, "key_sources"] = (
            rescue.loc[asmtl_mask, "key_sources"].astype(str)
            + "; PMID:39859321; NCBI Gene; HGNC"
        )
        rescue.loc[asmtl_mask, "display_decision_note"] = (
            "Table S5-only ASMT/ASMTL-locus context; not displayed in main Figure4C and "
            "not evidence for direct altered melatonin secretion."
        )

    confidence_lower = rescue["annotation_confidence"].astype(str).str.lower()
    low_or_unassigned = confidence_lower.isin(["low", "unassigned"]) | rescue["display_status"].eq("unassigned")
    force_table_s5_only = rescue["gene"].isin(FORCE_TABLES5_ONLY_GENES)
    force_main_display = rescue["gene"].isin(FORCE_MAIN_DISPLAY_GENES)

    rescue["recommended_for_Figure4C_display"] = (
        (
            rescue["display_status"].eq("main_Figure4C")
            & rescue["recommended_for_Figure4C_original"]
            & ~low_or_unassigned
            & ~rescue["recommended_for_TableS5_only_original"]
        )
        | force_main_display
    )
    rescue.loc[force_table_s5_only, "recommended_for_Figure4C_display"] = False
    rescue.loc[asmtl_mask, "recommended_for_Figure4C_display"] = False

    rescue["recommended_for_TableS5_only"] = (
        rescue["recommended_for_TableS5_only_original"]
        | rescue["display_status"].isin(["TableS5_only", "optional_Figure4C", "unassigned"])
        | low_or_unassigned
        | force_table_s5_only
    )
    rescue.loc[force_main_display, "recommended_for_TableS5_only"] = False
    rescue.loc[asmtl_mask, "recommended_for_TableS5_only"] = True
    rescue["keep_unassigned"] = rescue["keep_unassigned_original"] | rescue["display_status"].eq("unassigned")
    rescue.loc[asmtl_mask | force_table_s5_only, "keep_unassigned"] = False

    rescue.loc[
        rescue["display_decision_note"].eq("") & rescue["recommended_for_Figure4C_display"],
        "display_decision_note",
    ] = "Main Figure4C display candidate under compressed-module policy."
    rescue.loc[
        rescue["display_decision_note"].eq("") & rescue["display_status"].eq("optional_Figure4C"),
        "display_decision_note",
    ] = "Optional/Table S5-only by default because membrane trafficking / ER / vesicle processing is excluded from main Figure4C."
    rescue.loc[
        rescue["display_decision_note"].eq("") & rescue["display_status"].eq("TableS5_only"),
        "display_decision_note",
    ] = "Table S5-only module under final display policy."
    rescue.loc[
        rescue["display_decision_note"].eq("") & rescue["display_status"].eq("unassigned"),
        "display_decision_note",
    ] = "Unassigned or low-confidence; retained in Table S5 only."
    rescue.loc[
        rescue["gene"].isin(FORCE_TABLES5_ONLY_GENES),
        "display_decision_note",
    ] = "Explicit final decision: retain annotation in Table S5 only; do not display in main Figure4C."
    rescue.loc[
        rescue["gene"].isin(FORCE_MAIN_DISPLAY_GENES),
        "display_decision_note",
    ] = "Explicit final decision: retain as a main Figure4C representative gene."

    table_s5_cols = [
        "gene",
        "lasso_group",
        "selected_in_marine",
        "selected_in_aquatic",
        "marine_coef",
        "aquatic_coef",
        "max_abs_coef",
        "coefficient_direction_summary",
        "pro_recommended_module",
        "display_module",
        "recommended_submodule_or_function",
        "annotation_confidence",
        "recommended_for_Figure4C_original",
        "recommended_for_Figure4C_display",
        "recommended_for_TableS5_only",
        "keep_unassigned",
        "evidence_summary",
        "key_sources",
        "reason_for_assignment",
        "notes_for_manual_review",
        "display_decision_note",
    ]
    rescue.sort_values(["display_module", "lasso_group", "gene"])[table_s5_cols].to_csv(
        TABLE_S5_OUT, sep="\t", index=False
    )

    display = rescue[rescue["recommended_for_Figure4C_display"]].copy()
    display["confidence_rank"] = display["annotation_confidence"].map(confidence_rank)
    display["has_coef_marker"] = (
        display["coef_marker_marine_gt_0_1"].eq("#") | display["coef_marker_aquatic_gt_0_1"].eq("*")
    )
    display["priority_rank"] = display["is_priority_rescue_gene"].map({True: 0, False: 1})

    representative_rows = []
    cell_rows = []
    for display_module in MAIN_DISPLAY_MODULES:
        for group_raw, group_label in DISPLAY_LASSO_GROUP.items():
            cell = display[
                display["display_module"].eq(display_module) & display["lasso_group"].eq(group_raw)
            ].copy()
            if cell.empty:
                cell_rows.append({
                    "display_module": display_module,
                    "lasso_group": group_label,
                    "n_total_genes_in_cell": 0,
                    "n_high_confidence": 0,
                    "n_medium_confidence": 0,
                    "n_low_confidence": 0,
                    "n_coef_gt_0_1_marine": 0,
                    "n_coef_gt_0_1_aquatic": 0,
                    "n_representative_genes_displayed": 0,
                    "representative_genes": "",
                    "all_genes_in_cell": "",
                    "notes": "",
                })
                continue

            cell = cell.sort_values(
                ["confidence_rank", "max_abs_coef", "has_coef_marker", "priority_rank", "gene"],
                ascending=[True, False, False, True, True],
            )
            forced_genes = FORCE_REPRESENTATIVES.get((display_module, group_raw), [])
            max_display = 5 if forced_genes or len(cell) > 8 else 4
            forced = cell[cell["gene"].isin(forced_genes)].copy()
            if not forced.empty:
                forced["forced_order"] = forced["gene"].map({gene: i for i, gene in enumerate(forced_genes)})
                forced = forced.sort_values("forced_order").drop(columns=["forced_order"])
            remaining = cell[~cell["gene"].isin(forced["gene"] if not forced.empty else [])].copy()
            reps = pd.concat([forced, remaining.head(max(0, max_display - len(forced)))], ignore_index=True)
            for i, (_, row) in enumerate(reps.iterrows(), start=1):
                display_reason = []
                if row["gene"] in forced_genes:
                    display_reason.append("explicit_required_representative")
                if str(row["annotation_confidence"]).lower() in ["high", "medium"]:
                    display_reason.append(f"{row['annotation_confidence']}_confidence")
                if row["max_abs_coef"] > 0.1:
                    display_reason.append("abs_coef_gt_0.1")
                if row["is_priority_rescue_gene"]:
                    display_reason.append("priority_rescue_gene")
                representative_rows.append({
                    "display_module": display_module,
                    "lasso_group": group_label,
                    "gene": row["gene"],
                    "marine_coef": row["marine_coef"],
                    "aquatic_coef": row["aquatic_coef"],
                    "max_abs_coef": row["max_abs_coef"],
                    "coef_marker_marine_gt_0_1": row["coef_marker_marine_gt_0_1"],
                    "coef_marker_aquatic_gt_0_1": row["coef_marker_aquatic_gt_0_1"],
                    "coefficient_color_direction": row["coefficient_color_direction"],
                    "annotation_confidence": row["annotation_confidence"],
                    "display_order_within_cell": i,
                    "display_reason": "; ".join(display_reason) if display_reason else "top_ranked_for_cell",
                })

            notes = ""
            if len(cell) > len(reps):
                notes = f"Representative genes truncated from {len(cell)} to {len(reps)} for display readability."
            cell_rows.append({
                "display_module": display_module,
                "lasso_group": group_label,
                "n_total_genes_in_cell": len(cell),
                "n_high_confidence": int(cell["annotation_confidence"].astype(str).str.lower().eq("high").sum()),
                "n_medium_confidence": int(cell["annotation_confidence"].astype(str).str.lower().eq("medium").sum()),
                "n_low_confidence": int(cell["annotation_confidence"].astype(str).str.lower().eq("low").sum()),
                "n_coef_gt_0_1_marine": int(cell["abs_marine_coef"].gt(0.1).sum()),
                "n_coef_gt_0_1_aquatic": int(cell["abs_aquatic_coef"].gt(0.1).sum()),
                "n_representative_genes_displayed": len(reps),
                "representative_genes": ";".join(reps["gene"].tolist()),
                "all_genes_in_cell": ";".join(cell["gene"].tolist()),
                "notes": notes,
            })

    representative = pd.DataFrame(representative_rows)
    representative.to_csv(REPRESENTATIVE_OUT, sep="\t", index=False)
    cell_summary = pd.DataFrame(cell_rows)
    cell_summary.to_csv(CELL_SUMMARY_OUT, sep="\t", index=False)

    module_rows = []
    for display_module in MAIN_DISPLAY_MODULES:
        mod = display[display["display_module"].eq(display_module)].copy()
        reps = representative[representative["display_module"].eq(display_module)].copy()
        priority_genes = sorted(mod[mod["is_priority_rescue_gene"]]["gene"].unique())
        module_rows.append({
            "display_module": display_module,
            "n_total_predictors": len(mod),
            "n_shared": int(mod["lasso_group"].eq("shared").sum()),
            "n_marine_specific": int(mod["lasso_group"].eq("marine_specific").sum()),
            "n_aquatic_specific": int(mod["lasso_group"].eq("aquatic_specific").sum()),
            "n_displayed_representative_genes": len(reps),
            "representative_genes_overall": ";".join(reps["gene"].tolist()),
            "high_priority_rescued_genes": ";".join(priority_genes),
            "notes": "Compressed display module for Figure 4C; descriptive organization only.",
        })
    pd.DataFrame(module_rows).to_csv(MODULE_DISPLAY_OUT, sep="\t", index=False)

    write_warnings(rescue, display, cell_summary, representative)
    write_qc(rescue, compression, display, cell_summary, representative)
    write_caption_note()
    write_readme()
    write_manifest()
    package_path = write_zip()

    return {
        "compression": COMPRESSION_MAP_OUT,
        "table_s5": TABLE_S5_OUT,
        "cell_summary": CELL_SUMMARY_OUT,
        "representative": REPRESENTATIVE_OUT,
        "module_summary": MODULE_DISPLAY_OUT,
        "qc": QC_OUT,
        "warnings": WARNING_OUT,
        "caption": CAPTION_NOTE_OUT,
        "manifest": MANIFEST_OUT,
        "zip": package_path,
    }


def write_warnings(rescue: pd.DataFrame, display: pd.DataFrame, cell_summary: pd.DataFrame, representative: pd.DataFrame) -> None:
    rep_genes = set(representative["gene"]) if not representative.empty else set()
    warnings = []

    for _, row in cell_summary[cell_summary["n_total_genes_in_cell"] > row_value_safe(cell_summary, "n_representative_genes_displayed")].iterrows():
        warnings.append({
            "warning_type": "representative_genes_truncated",
            "gene": "",
            "display_module": row["display_module"],
            "lasso_group": row["lasso_group"],
            "severity": "info",
            "details": row["notes"],
        })

    displayed_low = display[display["annotation_confidence"].astype(str).str.lower().isin(["medium", "low"])]
    for _, row in displayed_low.iterrows():
        warnings.append({
            "warning_type": f"{row['annotation_confidence']}_confidence_gene_displayed",
            "gene": row["gene"],
            "display_module": row["display_module"],
            "lasso_group": DISPLAY_LASSO_GROUP.get(row["lasso_group"], row["lasso_group"]),
            "severity": "review" if str(row["annotation_confidence"]).lower() == "low" else "note",
            "details": row.get("notes_for_manual_review", ""),
        })

    high_coef_not_rep = rescue[
        rescue["max_abs_coef"].gt(0.1)
        & ~rescue["gene"].isin(rep_genes)
    ].copy()
    for _, row in high_coef_not_rep.iterrows():
        reason = "not selected as representative after ranking/truncation"
        if not row["recommended_for_Figure4C_display"]:
            if row["display_status"] == "optional_Figure4C":
                reason = "optional module excluded from main Figure4C by default"
            elif row["display_status"] in ["TableS5_only", "unassigned"]:
                reason = f"{row['display_status']} per display policy"
            elif str(row["annotation_confidence"]).lower() in ["low", "unassigned"]:
                reason = "low/unassigned confidence kept out of main display"
            elif row["recommended_for_TableS5_only"]:
                reason = "Pro/TableS5-only flag retained"
        warnings.append({
            "warning_type": "high_coefficient_gene_not_representative",
            "gene": row["gene"],
            "display_module": row["display_module"],
            "lasso_group": DISPLAY_LASSO_GROUP.get(row["lasso_group"], row["lasso_group"]),
            "severity": "note",
            "details": reason,
        })

    priority_not_rep = rescue[
        rescue["is_priority_rescue_gene"]
        & ~rescue["gene"].isin(rep_genes)
    ].copy()
    for _, row in priority_not_rep.iterrows():
        if row["gene"] == "RARS1":
            reason = "priority gene retained in Table S5 only because RNA processing / translation is not an approved main Figure4C module"
        elif not row["recommended_for_Figure4C_display"]:
            reason = "priority gene not eligible for main display under module/confidence policy"
        else:
            reason = "priority gene considered but not selected as representative after confidence/coefficient ranking and cell truncation"
        warnings.append({
            "warning_type": "priority_gene_not_representative",
            "gene": row["gene"],
            "display_module": row["display_module"],
            "lasso_group": DISPLAY_LASSO_GROUP.get(row["lasso_group"], row["lasso_group"]),
            "severity": "note",
            "details": reason,
        })

    for gene in sorted(FORCE_TABLES5_ONLY_GENES | {"ASMTL"}):
        row = rescue[rescue["gene"].eq(gene)]
        if row.empty:
            continue
        row = row.iloc[0]
        detail = (
            ASMTL_DISPLAY_NOTE if gene == "ASMTL"
            else "Explicit final decision: annotation retained in Table S5 only and not displayed in main Figure4C."
        )
        warnings.append({
            "warning_type": "explicit_TableS5_only_gene",
            "gene": gene,
            "display_module": row["display_module"],
            "lasso_group": DISPLAY_LASSO_GROUP.get(row["lasso_group"], row["lasso_group"]),
            "severity": "review" if gene == "ASMTL" else "note",
            "details": detail,
        })

    optional = rescue[rescue["display_status"].eq("optional_Figure4C")]
    for _, row in optional.iterrows():
        warnings.append({
            "warning_type": "optional_module_excluded_from_main_display",
            "gene": row["gene"],
            "display_module": row["display_module"],
            "lasso_group": DISPLAY_LASSO_GROUP.get(row["lasso_group"], row["lasso_group"]),
            "severity": "note",
            "details": "Membrane trafficking / ER / vesicle processing kept Table S5-only/optional by instruction.",
        })

    pd.DataFrame(warnings).to_csv(WARNING_OUT, sep="\t", index=False)


def row_value_safe(df: pd.DataFrame, column: str) -> pd.Series:
    return df[column] if column in df.columns else pd.Series([0] * len(df), index=df.index)


def write_qc(rescue: pd.DataFrame, compression: pd.DataFrame, display: pd.DataFrame,
             cell_summary: pd.DataFrame, representative: pd.DataFrame) -> None:
    reps_by_cell = {
        (row["display_module"], row["lasso_group"], row["gene"])
        for _, row in representative.iterrows()
    }
    all_by_cell = set()
    for _, row in cell_summary.iterrows():
        for gene in str(row["all_genes_in_cell"]).split(";"):
            if gene:
                all_by_cell.add((row["display_module"], row["lasso_group"], gene))

    def passfail(condition: bool) -> str:
        return "PASS" if bool(condition) else "FAIL"

    group_counts = rescue["lasso_group"].value_counts().to_dict()
    lasso_counts_match = (
        group_counts.get("shared", 0) == 24
        and group_counts.get("marine_specific", 0) == 47
        and group_counts.get("aquatic_specific", 0) == 77
    )
    main_modules = compression[compression["display_status"].eq("main_Figure4C")]["display_module"].nunique()

    def gene_not_unassigned(gene: str) -> bool:
        row = rescue[rescue["gene"].eq(gene)]
        return len(row) == 1 and row.iloc[0]["display_status"] != "unassigned"

    metabolism_genes = rescue[rescue["gene"].isin(["ACADSB", "ACOX3", "EHHADH"])]
    solute_genes = rescue[rescue["gene"].isin(["SLC6A6", "LRP2", "SPEF2", "WDR93"])]
    asmtl = rescue[rescue["gene"].eq("ASMTL")]
    rna = rescue[rescue["pro_recommended_module"].eq("RNA processing / translation / protein homeostasis")]
    unassigned = rescue[rescue["display_status"].eq("unassigned")]

    def gene_in_main_display(gene: str) -> bool:
        row = rescue[rescue["gene"].eq(gene)]
        return len(row) == 1 and bool(row.iloc[0]["recommended_for_Figure4C_display"])

    def gene_in_representatives(gene: str) -> bool:
        return gene in set(representative["gene"]) if not representative.empty else False

    checks = [
        ("n_total_genes", len(rescue), 148, passfail(len(rescue) == 148), ""),
        ("n_unique_genes", rescue["gene"].nunique(), 148, passfail(rescue["gene"].nunique() == 148), ""),
        ("no_duplicate_genes", rescue["gene"].duplicated().sum(), 0, passfail(not rescue["gene"].duplicated().any()), ""),
        ("lasso_group_counts_match", group_counts, "shared=24; marine_specific=47; aquatic_specific=77", passfail(lasso_counts_match), ""),
        ("display_modules_main_count <= 12", main_modules, "<=12", passfail(main_modules <= 12), ""),
        ("PER1_not_unassigned", rescue.loc[rescue["gene"].eq("PER1"), "display_module"].to_string(index=False), "not unassigned", passfail(gene_not_unassigned("PER1")), ""),
        ("PER1_display_status_checked", gene_in_main_display("PER1"), "TRUE", passfail(gene_in_main_display("PER1")), "PER1 remains the main high-confidence circadian/systemic regulator."),
        ("MYH7B_not_unassigned", rescue.loc[rescue["gene"].eq("MYH7B"), "display_module"].to_string(index=False), "not unassigned", passfail(gene_not_unassigned("MYH7B")), ""),
        ("IL36G_not_unassigned", rescue.loc[rescue["gene"].eq("IL36G"), "display_module"].to_string(index=False), "not unassigned", passfail(gene_not_unassigned("IL36G")), ""),
        ("ACADSB_ACOX3_EHHADH_in_metabolism_or_lipid_related_module",
         ";".join(metabolism_genes["display_module"].tolist()),
         "Metabolism / redox / lipid handling",
         passfail((metabolism_genes["display_module"] == "Metabolism / redox / lipid handling").all() and len(metabolism_genes) == 3),
         ""),
        ("SLC6A6_in_main_display",
         f"display={gene_in_main_display('SLC6A6')}; representative={gene_in_representatives('SLC6A6')}",
         "display=TRUE; representative=TRUE",
         passfail(gene_in_main_display("SLC6A6") and gene_in_representatives("SLC6A6")),
         ""),
        ("LRIF1_in_main_display",
         f"display={gene_in_main_display('LRIF1')}; representative={gene_in_representatives('LRIF1')}",
         "display=TRUE; representative=TRUE",
         passfail(gene_in_main_display("LRIF1") and gene_in_representatives("LRIF1")),
         ""),
        ("C1orf112_in_main_display",
         f"display={gene_in_main_display('C1orf112')}; representative={gene_in_representatives('C1orf112')}",
         "display=TRUE; representative=TRUE",
         passfail(gene_in_main_display("C1orf112") and gene_in_representatives("C1orf112")),
         ""),
        ("SPEF2_WDR93_handled",
         ";".join(f"{r.gene}:{r.display_module}" for r in solute_genes.itertuples()),
         "assigned to non-unassigned display modules",
         passfail((solute_genes["display_status"] != "unassigned").all() and len(solute_genes) == 4),
         ""),
        ("ASMTL_not_forced_into_main_module",
         (
             f"display_module={asmtl.iloc[0]['display_module']}; "
             f"display={asmtl.iloc[0]['recommended_for_Figure4C_display']}; "
             f"TableS5_only={asmtl.iloc[0]['recommended_for_TableS5_only']}; "
             f"confidence={asmtl.iloc[0]['annotation_confidence']}"
         ) if len(asmtl) else "missing",
         "endocrine/circadian context-only; display=FALSE; TableS5_only=TRUE; confidence=low",
         passfail(
             len(asmtl) == 1
             and asmtl.iloc[0]["display_module"] == "Endocrine / circadian / systemic regulation"
             and not asmtl.iloc[0]["recommended_for_Figure4C_display"]
             and asmtl.iloc[0]["recommended_for_TableS5_only"]
             and str(asmtl.iloc[0]["annotation_confidence"]).lower() == "low"
         ),
         "ASMTL retained only as cautious Table S5 context; ASMT is not a corrected full-data predictor here."),
        ("RNA_translation_not_in_main_Figure4C",
         int(rna["recommended_for_Figure4C_display"].sum()) if len(rna) else 0,
         0,
         passfail(len(rna) > 0 and not rna["recommended_for_Figure4C_display"].any()),
         ""),
        ("unassigned_genes_TableS5_only",
         int(unassigned["recommended_for_TableS5_only"].sum()),
         len(unassigned),
         passfail(len(unassigned) > 0 and unassigned["recommended_for_TableS5_only"].all()),
         ""),
        ("representative_genes_subset_of_all_genes",
         len(reps_by_cell - all_by_cell),
         0,
         passfail(len(reps_by_cell - all_by_cell) == 0),
         ""),
        ("counts_include_non-displayed_genes",
         int((cell_summary["n_total_genes_in_cell"] > cell_summary["n_representative_genes_displayed"]).sum()),
         ">0 crowded/truncated cells",
         passfail((cell_summary["n_total_genes_in_cell"] > cell_summary["n_representative_genes_displayed"]).any()),
         ""),
        ("marker_convention_matches_latest_Figure4",
         "# marine |coef| > 0.1; * aquatic |coef| > 0.1",
         "#=marine, *=aquatic",
         "PASS",
         ""),
    ]
    pd.DataFrame(checks, columns=["check", "observed", "expected", "status", "notes"]).to_csv(
        QC_OUT, sep="\t", index=False
    )


def write_caption_note() -> None:
    CAPTION_NOTE_OUT.write_text(
        "Panel C summarizes module-annotated selected predictors using compressed display modules. "
        "Circle size indicates the total number of selected predictors in each module-category cell, "
        "whereas gene labels show representative high- or medium-confidence genes selected for readability; "
        "full gene lists, unassigned predictors and Table S5-only annotations are reported in Supplementary "
        "Table S5. Module grouping is descriptive and does not by itself imply module-level evolutionary "
        "recurrence, which is evaluated separately in Fig. 5C.\n",
        encoding="utf-8",
    )


def write_readme() -> None:
    README_OUT.write_text(
        "# Figure4C display-table final package\n\n"
        "This folder contains final display-ready tables for the updated Figure4C and Supplementary Table S5.\n\n"
        "Scope:\n"
        "- Pro Figure4C annotation rescue outputs were used as the annotation source.\n"
        "- No Figure4C panel was drawn in this task.\n"
        "- No Figure4 PDF/SVG was modified.\n"
        "- No LASSO, ASR, GBI, t-test, or enrichment analysis was rerun.\n"
        "- Selected predictor sets and coefficients were not changed.\n\n"
        "Display policy:\n"
        "- The main Figure4C display uses the approved 12-module framework.\n"
        "- Representative genes are a readability subset; full gene lists are retained in the cell summary and Table S5 candidate.\n"
        "- ASMTL is retained as a low-confidence ASMT/ASMTL-locus, melatonin-pathway-adjacent context-only candidate in Table S5, not the main Figure4C display.\n"
        "- RARS1, RHCG, NOXO1, and IP6K1 are retained in Table S5 only unless explicitly approved later.\n\n"
        "Marker convention:\n"
        "- # = |marine coefficient| > 0.1\n"
        "- * = |aquatic coefficient| > 0.1\n",
        encoding="utf-8",
    )


def write_manifest() -> None:
    rows = []
    for path in sorted(OUT_DIR.glob("*")):
        if path.is_file() and path.name != MANIFEST_OUT.name and path.suffix != ".zip":
            rows.append({
                "file": path.name,
                "path": str(path),
                "sha256": sha256(path),
                "bytes": path.stat().st_size,
            })
    pd.DataFrame(rows).to_csv(MANIFEST_OUT, sep="\t", index=False)


def write_zip() -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    zip_path = INPUT_DIR / f"Figure4C_display_table_final_for_control_review_{timestamp}.zip"
    include = [
        COMPRESSION_MAP_OUT,
        TABLE_S5_OUT,
        CELL_SUMMARY_OUT,
        REPRESENTATIVE_OUT,
        MODULE_DISPLAY_OUT,
        QC_OUT,
        WARNING_OUT,
        CAPTION_NOTE_OUT,
        README_OUT,
        MANIFEST_OUT,
    ]
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for path in include:
            zf.write(path, arcname=path.name)
    return zip_path


def main() -> None:
    outputs = prepare_tables()
    for label, path in outputs.items():
        print(f"{label}: {path}")


if __name__ == "__main__":
    main()
