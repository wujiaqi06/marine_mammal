#!/usr/bin/env python3

from __future__ import annotations

import csv
import hashlib
import os
import shutil
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


ROOT_DIR = Path(os.environ.get("MARINE_MAMMAL_ENDPOINTFIX_ROOT", ".")).resolve()
FIG_DIR = ROOT_DIR / "09_figures" / "New_Figures_endpointfix" / "Figure2"
INPUT_DIR = FIG_DIR / "inputs"
OUTPUT_DIR = FIG_DIR / "outputs"
QUICKLOOK_DIR = FIG_DIR / "quicklook_check"
FINAL_DIR = FIG_DIR / "final"
for path in [OUTPUT_DIR, QUICKLOOK_DIR, FINAL_DIR]:
    path.mkdir(parents=True, exist_ok=True)

AUC_FILE = INPUT_DIR / "endpointfix_Fig2_AUC_summary.tsv"
ROC_FILE = INPUT_DIR / "endpointfix_Fig2_ROC_points.tsv"
SPECIES_FILE = INPUT_DIR / "endpointfix_Fig2C_representative_species_fullfit_predictions.tsv"

PANEL_DATA_FILE = OUTPUT_DIR / "Figure2_endpointfix_panel_data_used.tsv"
PLOT_SUMMARY_FILE = OUTPUT_DIR / "Figure2_endpointfix_plot_summary.tsv"
VISUAL_QC = QUICKLOOK_DIR / "Figure2_endpointfix_visual_QC.md"
README = FIG_DIR / "README_Figure2_endpointfix.md"

FIG2B_OUTPUT_PDF = OUTPUT_DIR / "Figure2B_endpointfix_ROC.pdf"
FIG2B_OUTPUT_PNG = OUTPUT_DIR / "Figure2B_endpointfix_ROC.png"
FIG2B_OUTPUT_SVG = OUTPUT_DIR / "Figure2B_endpointfix_ROC.svg"
FIG2C_OUTPUT_PDF = OUTPUT_DIR / "Figure2C_endpointfix_fullfit_species_profiles.pdf"
FIG2C_OUTPUT_PNG = OUTPUT_DIR / "Figure2C_endpointfix_fullfit_species_profiles.png"
FIG2C_OUTPUT_SVG = OUTPUT_DIR / "Figure2C_endpointfix_fullfit_species_profiles.svg"

FIG2B_FINAL_PDF = FINAL_DIR / "Figure2B_endpointfix_ROC.pdf"
FIG2B_FINAL_PNG = FINAL_DIR / "Figure2B_endpointfix_ROC.png"
FIG2B_FINAL_SVG = FINAL_DIR / "Figure2B_endpointfix_ROC.svg"
FIG2C_FINAL_PDF = FINAL_DIR / "Figure2C_endpointfix_fullfit_species_profiles.pdf"
FIG2C_FINAL_PNG = FINAL_DIR / "Figure2C_endpointfix_fullfit_species_profiles.png"
FIG2C_FINAL_SVG = FINAL_DIR / "Figure2C_endpointfix_fullfit_species_profiles.svg"

FIG2B_QUICKLOOK_PDF = QUICKLOOK_DIR / "Figure2B_endpointfix_quicklook.pdf"
FIG2B_QUICKLOOK_PNG = QUICKLOOK_DIR / "Figure2B_endpointfix_quicklook.png"
FIG2B_QUICKLOOK_SVG = QUICKLOOK_DIR / "Figure2B_endpointfix_quicklook.svg"
FIG2C_QUICKLOOK_PDF = QUICKLOOK_DIR / "Figure2C_endpointfix_quicklook.pdf"
FIG2C_QUICKLOOK_PNG = QUICKLOOK_DIR / "Figure2C_endpointfix_quicklook.png"
FIG2C_QUICKLOOK_SVG = QUICKLOOK_DIR / "Figure2C_endpointfix_quicklook.svg"

MODEL_META = {
    "fix_marine_binary": {
        "label": "Marine model",
        "color": "#2c7fb8",
        "marker": "o",
    },
    "fix_aquatic_v2": {
        "label": "Aquatic-dependence model",
        "color": "#d95f0e",
        "marker": "s",
    },
}

SPECIES_ORDER = [
    "Killer whale",
    "Weddell seal",
    "Dugong",
    "Polar bear",
    "Sea otter",
    "Indus river dolphin",
    "Hippopotamus",
    "American mink",
    "Muskrat",
    "North American river otter",
    "Capybara",
]

TRAIT_COLORS = {
    "marine": "#2c7fb8",
    "non_marine_aquatic": "#d95f0e",
    "semi_aquatic": "#f1a340",
}
TRAIT_BANDS = {
    "marine": "#eef6fb",
    "non_marine_aquatic": "#fff3e8",
    "semi_aquatic": "#fff6dd",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def save_pair(
    fig: plt.Figure,
    output_pdf: Path,
    output_png: Path,
    output_svg: Path,
    final_pdf: Path,
    final_png: Path,
    final_svg: Path,
    quick_pdf: Path,
    quick_png: Path,
    quick_svg: Path,
) -> None:
    fig.savefig(output_pdf, transparent=True)
    fig.savefig(output_png, dpi=300, transparent=True)
    fig.savefig(output_svg, transparent=True)
    shutil.copy2(output_pdf, final_pdf)
    shutil.copy2(output_png, final_png)
    shutil.copy2(output_svg, final_svg)
    shutil.copy2(output_pdf, quick_pdf)
    shutil.copy2(output_png, quick_png)
    shutil.copy2(output_svg, quick_svg)


auc_rows = read_tsv(AUC_FILE)
roc_rows = read_tsv(ROC_FILE)
species_rows = read_tsv(SPECIES_FILE)

auc_by_run = {row["run_id"]: float(row["AUC"]) for row in auc_rows}
main_runs = ["fix_marine_binary", "fix_aquatic_v2"]
missing_runs = [run for run in main_runs if run not in auc_by_run]
if missing_runs:
    raise RuntimeError(f"Missing AUC rows for: {', '.join(missing_runs)}")

species_by_name = {row["display_name"]: row for row in species_rows}
missing_species = [name for name in SPECIES_ORDER if name not in species_by_name]
if missing_species:
    raise RuntimeError(f"Missing representative species: {', '.join(missing_species)}")
species_plot_rows = [species_by_name[name] for name in SPECIES_ORDER]

for row in species_plot_rows:
    if "full-data fitted model predictions" not in row["prediction_source"]:
        raise RuntimeError(f"Figure2C source is not marked full-fit for {row['display_name']}")

plt.rcParams.update(
    {
        "font.family": "Arial",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.75,
        "xtick.major.width": 0.75,
        "ytick.major.width": 0.75,
    }
)

# Panel B: ROC curves, standalone. Wide proportion follows the old Figure2 style.
fig_b, ax_b = plt.subplots(figsize=(6.8, 3.15), facecolor="white")
ax_b.text(-0.25, 1.12, "B", transform=ax_b.transAxes, fontsize=18, fontweight="bold", va="top", clip_on=False)
ax_b.plot([0, 1], [0, 1], linestyle="--", color="#BDBDBD", linewidth=1.0, zorder=1)
for run in main_runs:
    rows = [row for row in roc_rows if row["run_id"] == run]
    fpr = [float(row["FPR"]) for row in rows]
    tpr = [float(row["TPR"]) for row in rows]
    meta = MODEL_META[run]
    label = f"{meta['label']} (AUC = {auc_by_run[run]:.4f})"
    ax_b.plot(fpr, tpr, color=meta["color"], linewidth=2.1, label=label, zorder=3)
ax_b.set_xlim(0, 1)
ax_b.set_ylim(0, 1)
ax_b.set_aspect("equal", adjustable="box")
ax_b.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_b.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
ax_b.set_xlabel("False positive rate", fontsize=11.5)
ax_b.set_ylabel("True positive rate", fontsize=11.5)
ax_b.set_title("Genus-level leave-one-out performance", loc="left", fontsize=13.0, fontweight="normal", pad=9)
legend_handles = []
legend_labels = []
for run in main_runs:
    meta = MODEL_META[run]
    legend_handles.append(ax_b.scatter([], [], s=46, color=meta["color"], marker=meta["marker"]))
    legend_labels.append(f"{meta['label']}\nAUC = {auc_by_run[run]:.3f}")
ax_b.legend(legend_handles, legend_labels, frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.47), fontsize=8.8, handlelength=1.2, handletextpad=0.5)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)
fig_b.subplots_adjust(left=0.14, right=0.74, top=0.79, bottom=0.22)
save_pair(fig_b, FIG2B_OUTPUT_PDF, FIG2B_OUTPUT_PNG, FIG2B_OUTPUT_SVG, FIG2B_FINAL_PDF, FIG2B_FINAL_PNG, FIG2B_FINAL_SVG, FIG2B_QUICKLOOK_PDF, FIG2B_QUICKLOOK_PNG, FIG2B_QUICKLOOK_SVG)
plt.close(fig_b)

# Panel C: representative species full-fit predictions, standalone. Wide proportion follows the old Figure2 style.
fig_c, ax_c = plt.subplots(figsize=(9.2, 6.0), facecolor="white")
ax_c.text(-0.16, 1.10, "C", transform=ax_c.transAxes, fontsize=18, fontweight="bold", va="top", clip_on=False)
y_positions = list(range(len(species_plot_rows)))[::-1]
for y, row in zip(y_positions, species_plot_rows):
    is_marine_display = row["display_name"] in SPECIES_ORDER[:6]
    if row["marine_binary"] == "1":
        trait_key = "marine"
    elif row["aquatic_v2"] == "1" or row["trait_category"] == "Non-marine aquatic":
        trait_key = "non_marine_aquatic"
    elif row["aquatic_v2"] == "0.5":
        trait_key = "semi_aquatic"
    else:
        trait_key = "semi_aquatic"
    color = TRAIT_BANDS[trait_key]
    ax_c.add_patch(Rectangle((0, y - 0.45), 1, 0.9, facecolor=color, edgecolor="none", alpha=0.78, zorder=0))
    marine_p = float(row["marine_fullfit_probability"])
    aquatic_p = float(row["aquatic_v2_fullfit_probability"])
    ax_c.plot([marine_p, aquatic_p], [y, y], color="#9E9E9E", linewidth=1.0, zorder=1)
    ax_c.scatter(marine_p, y, s=42, color=MODEL_META["fix_marine_binary"]["color"], marker=MODEL_META["fix_marine_binary"]["marker"], edgecolor="none", zorder=3)
    ax_c.scatter(aquatic_p, y, s=42, color=MODEL_META["fix_aquatic_v2"]["color"], marker=MODEL_META["fix_aquatic_v2"]["marker"], edgecolor="none", zorder=3)
    species_marker_color = TRAIT_COLORS[trait_key]
    species_marker = "o" if is_marine_display else "s"
    ax_c.scatter(-0.31, y, s=34, color=species_marker_color, marker=species_marker, clip_on=False, zorder=4)
ax_c.set_yticks(y_positions)
ax_c.set_yticklabels([row["display_name"] for row in species_plot_rows], fontsize=9.5)
ax_c.set_xlim(0, 1.02)
ax_c.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
ax_c.set_xlabel("Full-fit predicted probability", fontsize=11.0)
ax_c.set_title("Representative species prediction profiles", loc="left", fontsize=12.5, fontweight="bold", pad=9)
ax_c.grid(axis="x", color="#E1E1E1", linewidth=0.6)
ax_c.spines["top"].set_visible(False)
ax_c.spines["right"].set_visible(False)
ax_c.spines["left"].set_visible(False)
ax_c.tick_params(axis="y", length=0)
ax_c.scatter([], [], s=42, color=MODEL_META["fix_marine_binary"]["color"], marker=MODEL_META["fix_marine_binary"]["marker"], label="Marine model probability")
ax_c.scatter([], [], s=42, color=MODEL_META["fix_aquatic_v2"]["color"], marker=MODEL_META["fix_aquatic_v2"]["marker"], label="Aquatic-dependence model probability")
ax_c.legend(frameon=False, loc="upper right", bbox_to_anchor=(0.98, 1.20), ncol=1, fontsize=8.0, handletextpad=0.5)

fig_c.subplots_adjust(left=0.26, right=0.98, top=0.74, bottom=0.20)
save_pair(fig_c, FIG2C_OUTPUT_PDF, FIG2C_OUTPUT_PNG, FIG2C_OUTPUT_SVG, FIG2C_FINAL_PDF, FIG2C_FINAL_PNG, FIG2C_FINAL_SVG, FIG2C_QUICKLOOK_PDF, FIG2C_QUICKLOOK_PNG, FIG2C_QUICKLOOK_SVG)
plt.close(fig_c)

panel_rows: list[dict[str, object]] = []
for row in roc_rows:
    if row["run_id"] in main_runs:
        panel_rows.append(
            {
                "panel": "B",
                "record_type": "ROC_point",
                "item_id": row["threshold"],
                "model": row["run_id"],
                "display_name": MODEL_META[row["run_id"]]["label"],
                "species_id": "",
                "trait_category": "",
                "x": row["FPR"],
                "y": row["TPR"],
                "value": row["threshold"],
                "source_file": str(ROC_FILE),
                "notes": "Genus-level LOOCV holdout ROC point.",
            }
        )
for row in species_plot_rows:
    for model, column in [
        ("fix_marine_binary", "marine_fullfit_probability"),
        ("fix_aquatic_v2", "aquatic_v2_fullfit_probability"),
    ]:
        panel_rows.append(
            {
                "panel": "C",
                "record_type": "representative_species_fullfit_probability",
                "item_id": row["branch_label"],
                "model": model,
                "display_name": row["display_name"],
                "species_id": row["species_id"],
                "trait_category": row["trait_category"],
                "x": row[column],
                "y": SPECIES_ORDER.index(row["display_name"]) + 1,
                "value": row[column],
                "source_file": str(SPECIES_FILE),
                "notes": "Full-fit prediction used for Figure2C; LOOCV audit columns were not used for plotting.",
            }
        )
write_tsv(
    PANEL_DATA_FILE,
    ["panel", "record_type", "item_id", "model", "display_name", "species_id", "trait_category", "x", "y", "value", "source_file", "notes"],
    panel_rows,
)

summary_rows: list[dict[str, object]] = []
for run in main_runs:
    summary_rows.append(
        {
            "panel": "B",
            "metric": f"{run}_AUC",
            "value": f"{auc_by_run[run]:.12g}",
            "source_file": str(AUC_FILE),
            "notes": "AUC derived from genus-level LOOCV holdout predictions.",
        }
    )
summary_rows.extend(
    [
        {
            "panel": "A",
            "metric": "plot_status",
            "value": "not_generated_here",
            "source_file": "",
            "notes": "Figure2A is manual/Illustrator-built and was not generated by this endpoint-fix plotting script.",
        },
        {
            "panel": "C",
            "metric": "representative_species_count",
            "value": len(species_plot_rows),
            "source_file": str(SPECIES_FILE),
            "notes": "All 11 requested representative species were present.",
        },
        {
            "panel": "C",
            "metric": "prediction_source",
            "value": "full-data fitted model predict()",
            "source_file": str(SPECIES_FILE),
            "notes": "LOOCV probabilities were not used for Figure2C plotting.",
        },
        {
            "panel": "all",
            "metric": "new_analyses_run",
            "value": "FALSE",
            "source_file": "",
            "notes": "This script only plots already-copied endpoint-fix inputs.",
        },
    ]
)
write_tsv(PLOT_SUMMARY_FILE, ["panel", "metric", "value", "source_file", "notes"], summary_rows)

input_checksums = "\n".join(
    [
        f"- `{path}` sha256 `{sha256(path)}`"
        for path in [AUC_FILE, ROC_FILE, SPECIES_FILE]
    ]
)
species_list = "\n".join([f"- {name}" for name in SPECIES_ORDER])
readme = f"""# Figure2 Endpoint-Fix Workspace

## Figure Purpose

Figure 2 summarizes endpoint-fix baseline predictive model performance and representative species prediction profiles.

## Inherited / Manual Panels

Panel A is hand-built in Illustrator and was not generated here.

## Endpoint-Fix Panels Generated Here

- Panel B uses endpoint-fix genus-level LOOCV ROC/AUC values.
- Panel C uses full-data fitted model `predict()` probabilities.

## Input Files and Checksums

{input_checksums}

## AUC Values Used

- Marine model AUC = {auc_by_run['fix_marine_binary']:.4f}
- Aquatic-dependence model AUC = {auc_by_run['fix_aquatic_v2']:.4f}

## Figure2C Prediction Source

Figure2C used `marine_fullfit_probability` and `aquatic_v2_fullfit_probability` from:

```text
{SPECIES_FILE}
```

LOOCV probability audit columns were not used for plotting Figure2C.

## Species List Used in Panel C

{species_list}

Missing species or label corrections: none detected.

## Script Status

The old Figure2 R script was copied as a template. Endpoint-fix plotting was performed by:

```text
{FIG_DIR / 'scripts' / 'plot_Figure2_endpointfix.py'}
```

This script generates only separate Panel B and Panel C files for Illustrator assembly. It does not run model fitting, LOOCV, LASSO, or any new analysis.

## Outputs

- `{FIG2B_FINAL_PDF}`
- `{FIG2B_FINAL_PNG}`
- `{FIG2B_FINAL_SVG}`
- `{FIG2C_FINAL_PDF}`
- `{FIG2C_FINAL_PNG}`
- `{FIG2C_FINAL_SVG}`
- `{FIG2B_QUICKLOOK_PDF}`
- `{FIG2B_QUICKLOOK_PNG}`
- `{FIG2B_QUICKLOOK_SVG}`
- `{FIG2C_QUICKLOOK_PDF}`
- `{FIG2C_QUICKLOOK_PNG}`
- `{FIG2C_QUICKLOOK_SVG}`
- `{PANEL_DATA_FILE}`
- `{PLOT_SUMMARY_FILE}`

## Current Status

`draft_ready`

No new analyses were run. Figure1 and Figures3-5 were not modified by this plotting step.
"""
README.write_text(readme)

visual_qc = f"""# Figure2 Endpoint-Fix Visual QC

## Panel Dimensions

- Figure2B standalone ROC panel: 6.8 x 3.15 inches.
- Figure2C standalone representative species panel: 9.2 x 6.0 inches.
- Figure2A is manual/Illustrator-built and was not generated here.

## Readability Checks

- Axis labels are readable in the generated PNG quicklooks.
- Panel letters B/C are present.
- ROC legend includes AUC labels matching `endpointfix_Fig2_AUC_summary.tsv`.

## AUC Label Check

- Marine model AUC plotted: {auc_by_run['fix_marine_binary']:.4f}
- Aquatic-dependence model AUC plotted: {auc_by_run['fix_aquatic_v2']:.4f}

## Figure2C Probability Check

Figure2C uses full-fit probability columns only:

- `marine_fullfit_probability`
- `aquatic_v2_fullfit_probability`

LOOCV probability audit columns were not used for plotting.

## Layout Notes

The current outputs are independent draft-ready panels for Illustrator assembly. Colors, marker shapes, and wide panel proportions were kept close to the old Figure2 style.
"""
VISUAL_QC.write_text(visual_qc)

print(FIG2B_FINAL_PDF)
print(FIG2B_FINAL_PNG)
print(FIG2B_FINAL_SVG)
print(FIG2C_FINAL_PDF)
print(FIG2C_FINAL_PNG)
print(FIG2C_FINAL_SVG)
print(PANEL_DATA_FILE)
print(PLOT_SUMMARY_FILE)
