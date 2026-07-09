#!/usr/bin/env python3
"""Assemble corrected Figure 4 panels into a review PDF/SVG.

The combined PNG is intentionally not generated because earlier review packages
hit stale-raster contamination. Individual panel PNGs remain available.
"""

from __future__ import annotations

import hashlib
import html
import os
from pathlib import Path
from urllib.parse import quote

import pandas as pd
from PIL import Image
from pypdf import PdfReader, PdfWriter, PageObject, Transformation
from reportlab.pdfgen import canvas
from reportlab.lib.colors import black


SCRIPT_DIR = Path(__file__).resolve().parent
FIGURE_DIR = SCRIPT_DIR.parent
FINAL_DIR = FIGURE_DIR / "final"
QC_DIR = FIGURE_DIR / "qc"
INPUT_DIR = FIGURE_DIR / "inputs"
FINAL_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

PANEL_A = FIGURE_DIR / "Figure4A_outputs" / "Figure4A_gene_set_partition.pdf"
PANEL_B = FIGURE_DIR / "Figure4B_outputs" / "Figure4B_coefficient_architecture.pdf"
PANEL_C = FIGURE_DIR / "Figure4C_outputs" / "Figure4C_module_annotated_partition.pdf"
PANEL_A_PNG = FIGURE_DIR / "Figure4A_outputs" / "Figure4A_gene_set_partition.png"
PANEL_B_PNG = FIGURE_DIR / "Figure4B_outputs" / "Figure4B_coefficient_architecture.png"
PANEL_C_PNG = FIGURE_DIR / "Figure4C_outputs" / "Figure4C_module_annotated_partition.png"
PANEL_A_SVG = FIGURE_DIR / "Figure4A_outputs" / "Figure4A_gene_set_partition.svg"
PANEL_B_SVG = FIGURE_DIR / "Figure4B_outputs" / "Figure4B_coefficient_architecture.svg"
PANEL_C_SVG = FIGURE_DIR / "Figure4C_outputs" / "Figure4C_module_annotated_partition.svg"

COMBINED_PDF = FINAL_DIR / "Figure4_corrected_full_data_old_design.pdf"
COMBINED_SVG = FINAL_DIR / "Figure4_corrected_full_data_old_design.svg"
CAPTION_NOTE = FINAL_DIR / "Figure4_caption_evidence_note.md"
QC_TABLE = QC_DIR / "Figure4_corrected_old_design_check.tsv"
README = FIGURE_DIR / "README_Figure4_corrected_old_design.md"

POINTS_PER_INCH = 72
PAGE_W = 15.8 * POINTS_PER_INCH
PAGE_H = 13.8 * POINTS_PER_INCH


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def merge_panel(blank: PageObject, panel_path: Path, x: float, y: float, w: float, h: float) -> None:
    panel = PdfReader(str(panel_path)).pages[0]
    panel_w = float(panel.mediabox.width)
    panel_h = float(panel.mediabox.height)
    scale = min(w / panel_w, h / panel_h)
    tx = x + (w - panel_w * scale) / 2
    ty = y + (h - panel_h * scale) / 2
    panel.add_transformation(Transformation().scale(scale).translate(tx, ty))
    blank.merge_page(panel)


def make_label_overlay() -> Path:
    overlay = FINAL_DIR / "_Figure4_panel_labels_overlay.pdf"
    c = canvas.Canvas(str(overlay), pagesize=(PAGE_W, PAGE_H))
    c.setFillColor(black)
    c.setFont("Helvetica-Bold", 18)
    c.drawString(0.25 * POINTS_PER_INCH, PAGE_H - 0.45 * POINTS_PER_INCH, "A")
    c.drawString(5.05 * POINTS_PER_INCH, PAGE_H - 0.45 * POINTS_PER_INCH, "B")
    c.drawString(0.25 * POINTS_PER_INCH, 7.15 * POINTS_PER_INCH, "C")
    c.save()
    return overlay


def assemble_pdf() -> None:
    blank = PageObject.create_blank_page(width=PAGE_W, height=PAGE_H)

    # Layout: A and B across the top, C across the bottom.
    merge_panel(
        blank,
        PANEL_A,
        x=0.35 * POINTS_PER_INCH,
        y=7.40 * POINTS_PER_INCH,
        w=4.45 * POINTS_PER_INCH,
        h=5.85 * POINTS_PER_INCH,
    )
    merge_panel(
        blank,
        PANEL_B,
        x=4.95 * POINTS_PER_INCH,
        y=7.30 * POINTS_PER_INCH,
        w=10.55 * POINTS_PER_INCH,
        h=5.95 * POINTS_PER_INCH,
    )
    merge_panel(
        blank,
        PANEL_C,
        x=0.35 * POINTS_PER_INCH,
        y=0.25 * POINTS_PER_INCH,
        w=15.10 * POINTS_PER_INCH,
        h=6.85 * POINTS_PER_INCH,
    )

    overlay = PdfReader(str(make_label_overlay())).pages[0]
    blank.merge_page(overlay)

    writer = PdfWriter()
    writer.add_page(blank)
    with COMBINED_PDF.open("wb") as f:
        writer.write(f)

    temp_overlay = FINAL_DIR / "_Figure4_panel_labels_overlay.pdf"
    if temp_overlay.exists():
        temp_overlay.unlink()


def svg_image_href(path: Path) -> str:
    rel = os.path.relpath(path, start=COMBINED_SVG.parent)
    return quote(str(rel))


def assemble_svg() -> None:
    width = 15.8 * 96
    height = 13.8 * 96
    pieces = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0f}" height="{height:.0f}" viewBox="0 0 {width:.0f} {height:.0f}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
        '<text x="24" y="40" font-family="Arial, Helvetica, sans-serif" font-size="24" font-weight="700">A</text>',
        '<text x="485" y="40" font-family="Arial, Helvetica, sans-serif" font-size="24" font-weight="700">B</text>',
        '<text x="24" y="642" font-family="Arial, Helvetica, sans-serif" font-size="24" font-weight="700">C</text>',
        f'<image href="{html.escape(svg_image_href(PANEL_A_SVG))}" x="34" y="54" width="427" height="562" preserveAspectRatio="xMidYMid meet"/>',
        f'<image href="{html.escape(svg_image_href(PANEL_B_SVG))}" x="475" y="54" width="1013" height="571" preserveAspectRatio="xMidYMid meet"/>',
        f'<image href="{html.escape(svg_image_href(PANEL_C_SVG))}" x="34" y="656" width="1450" height="658" preserveAspectRatio="xMidYMid meet"/>',
        "</svg>",
    ]
    COMBINED_SVG.write_text("\n".join(pieces) + "\n", encoding="utf-8")


def write_caption_note() -> None:
    CAPTION_NOTE.write_text(
        "Figure 4 summarizes the final full-data LASSO architecture after corrected preprocessing. "
        "It is not a cross-validated validation panel. Module grouping is a descriptive organization "
        "of module-annotated selected predictors and does not by itself imply module-level evolutionary "
        "recurrence.\n",
        encoding="utf-8",
    )


def write_readme() -> None:
    README.write_text(
        "# Figure 4 Corrected Full-Data Old-Design Outputs\n\n"
        "This folder rebuilds Figure 4 using the old visual design family while replacing all legacy "
        "numbers and gene sets with Phase 12A corrected-preprocessing full-data baseline LASSO outputs.\n\n"
        "Corrected values used: marine predictors = 71, binary aquatic-dependence predictors = 101, "
        "shared = 24, marine-specific = 47, aquatic-specific = 77, union total = 148.\n\n"
        "Panel C lists all module-annotated gene names in the plot table and rendered panel. "
        "Unannotated selected predictors are retained in the gene-level input table but excluded from "
        "module rows: shared = 3, marine-specific = 24, aquatic-specific = 45.\n\n"
        "The combined PNG is intentionally not generated because previous package review found stale "
        "raster contamination in combined Figure 4 PNG exports. Use the combined PDF/SVG and individual "
        "panel PNG/PDF/SVG files for review and Illustrator assembly.\n\n"
        "Caption/evidence note: Figure 4 summarizes the final full-data LASSO architecture after "
        "corrected preprocessing. It is not a cross-validated validation panel. Module grouping is a "
        "descriptive organization of module-annotated selected predictors and does not by itself imply "
        "module-level evolutionary recurrence.\n",
        encoding="utf-8",
    )


def get_png_dim(path: Path) -> str:
    im = Image.open(path)
    return f"{im.width}x{im.height}"


def write_qc() -> None:
    gene = pd.read_csv(INPUT_DIR / "Table1_gene_level_working_table.tsv", sep="\t")
    unannotated = pd.read_csv(INPUT_DIR / "Figure4_unannotated_predictor_counts.tsv", sep="\t")
    module_plot = pd.read_csv(FIGURE_DIR / "Figure4C_outputs" / "Figure4C_plot_table.tsv", sep="\t")

    marine_n = int(gene["selected_in_marine"].sum())
    aquatic_n = int(gene["selected_in_aquatic"].sum())
    shared_n = int(((gene["selected_in_marine"]) & (gene["selected_in_aquatic"])).sum())
    marine_specific_n = int(((gene["selected_in_marine"]) & (~gene["selected_in_aquatic"])).sum())
    aquatic_specific_n = int(((~gene["selected_in_marine"]) & (gene["selected_in_aquatic"])).sum())
    union_n = len(gene)

    unannotated_key = {
        "Shared predictors": "shared",
        "Marine-specific predictors": "marine_specific",
        "Aquatic-specific predictors": "aquatic_specific",
    }
    unannotated_counts = {
        unannotated_key.get(row["display_category"], row["display_category"]): int(row["gene_count"])
        for _, row in unannotated.iterrows()
    }

    checks = [
        ("marine_predictors", marine_n, 71, marine_n == 71),
        ("binary_aquatic_predictors", aquatic_n, 101, aquatic_n == 101),
        ("shared_predictors", shared_n, 24, shared_n == 24),
        ("marine_specific_predictors", marine_specific_n, 47, marine_specific_n == 47),
        ("aquatic_specific_predictors", aquatic_specific_n, 77, aquatic_specific_n == 77),
        ("union_total", union_n, 148, union_n == 148),
        ("shared_unannotated", unannotated_counts.get("shared", -1), 3, unannotated_counts.get("shared", -1) == 3),
        ("marine_specific_unannotated", unannotated_counts.get("marine_specific", -1), 24, unannotated_counts.get("marine_specific", -1) == 24),
        ("aquatic_specific_unannotated", unannotated_counts.get("aquatic_specific", -1), 45, unannotated_counts.get("aquatic_specific", -1) == 45),
        ("panel_A_pdf_exists", PANEL_A.exists(), True, PANEL_A.exists()),
        ("panel_A_svg_exists", PANEL_A_SVG.exists(), True, PANEL_A_SVG.exists()),
        ("panel_A_png_exists", PANEL_A_PNG.exists(), True, PANEL_A_PNG.exists()),
        ("panel_B_pdf_exists", PANEL_B.exists(), True, PANEL_B.exists()),
        ("panel_B_svg_exists", PANEL_B_SVG.exists(), True, PANEL_B_SVG.exists()),
        ("panel_B_png_exists", PANEL_B_PNG.exists(), True, PANEL_B_PNG.exists()),
        ("panel_C_pdf_exists", PANEL_C.exists(), True, PANEL_C.exists()),
        ("panel_C_svg_exists", PANEL_C_SVG.exists(), True, PANEL_C_SVG.exists()),
        ("panel_C_png_exists", PANEL_C_PNG.exists(), True, PANEL_C_PNG.exists()),
        ("combined_pdf_exists", COMBINED_PDF.exists(), True, COMBINED_PDF.exists()),
        ("combined_svg_exists", COMBINED_SVG.exists(), True, COMBINED_SVG.exists()),
        ("combined_png_intentionally_not_generated", not (FINAL_DIR / "Figure4_corrected_full_data_old_design.png").exists(), True, not (FINAL_DIR / "Figure4_corrected_full_data_old_design.png").exists()),
        ("panel_C_gene_names_all_listed", bool((module_plot.query("gene_count > 0")["gene_names"].astype(str).str.len() > 0).all()), True, bool((module_plot.query("gene_count > 0")["gene_names"].astype(str).str.len() > 0).all())),
    ]

    contaminant_terms = [
        "Inner Product",
        "Length of a Vector",
        "Distance in Rn",
        "Quality Control Engineer",
        "discharged effluent",
        "compliance standards",
    ]
    searchable_files = [
        COMBINED_SVG,
        README,
        CAPTION_NOTE,
        PLOT_TABLE_FILE := FIGURE_DIR / "Figure4C_outputs" / "Figure4C_plot_table.tsv",
    ]
    searchable_text = "\n".join(
        path.read_text(encoding="utf-8", errors="ignore")
        for path in searchable_files
        if path.exists()
    )
    contaminant_absent = not any(term in searchable_text for term in contaminant_terms)
    checks.append(("known_contaminant_strings_absent", contaminant_absent, True, contaminant_absent))

    rows = []
    for check, observed, expected, passed in checks:
        rows.append({
            "check": check,
            "observed": observed,
            "expected": expected,
            "status": "PASS" if passed else "FAIL",
            "notes": "",
        })

    for path in [PANEL_A_PNG, PANEL_B_PNG, PANEL_C_PNG]:
        rows.append({
            "check": f"{path.stem}_dimensions",
            "observed": get_png_dim(path),
            "expected": "image readable",
            "status": "PASS",
            "notes": f"sha256={sha256(path)}",
        })

    pd.DataFrame(rows).to_csv(QC_TABLE, sep="\t", index=False)


def main() -> None:
    for path in [PANEL_A, PANEL_B, PANEL_C, PANEL_A_SVG, PANEL_B_SVG, PANEL_C_SVG, PANEL_A_PNG, PANEL_B_PNG, PANEL_C_PNG]:
        if not path.exists():
            raise FileNotFoundError(path)
    assemble_pdf()
    assemble_svg()
    write_caption_note()
    write_readme()
    write_qc()
    print(f"Saved: {COMBINED_PDF}")
    print(f"Saved: {COMBINED_SVG}")
    print(f"Saved: {CAPTION_NOTE}")
    print(f"Saved: {QC_TABLE}")
    print(f"Saved: {README}")
    print("Combined PNG intentionally not generated.")


if __name__ == "__main__":
    main()
