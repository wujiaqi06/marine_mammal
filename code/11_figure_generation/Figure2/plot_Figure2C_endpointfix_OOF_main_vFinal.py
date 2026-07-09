#!/usr/bin/env python3
"""Draw standalone endpoint-fix Figure 2C from existing OOF predictions only.

This script only plots existing all-species genus-level LOOCV out-of-fold
predictions. It does not refit LASSO, rerun prediction, or edit other panels.
"""

from __future__ import annotations

import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


ROOT = Path(os.environ.get("MARINE_MAMMAL_ENDPOINTFIX_ROOT", ".")).resolve()
FIG2 = ROOT / "09_figures" / "New_Figures_endpointfix" / "Figure2"
INPUT = ROOT / "10_reviewer_risk_controls" / "07_pro_decision_integration" / "Fig2C_main_candidates" / "Figure2C_all_species_OOF_distribution_main_candidate_data.tsv"

OUTPUTS = FIG2 / "outputs"
FINAL = FIG2 / "final"
QUICK = FIG2 / "quicklook_check"
QC = ROOT / "10_reviewer_risk_controls" / "08_integrated_Result2.1_2.3_QC" / "figure2_outputs"

BLUE = "#2c7fb8"
ORANGE = "#d95f0e"
SEMI = "#f1a340"
GREY = "#bdbdbd"
TEXT_GREY = "#4d4d4d"
CATEGORY_ORDER = ["terrestrial", "semi-aquatic", "non-marine aquatic", "marine"]
COLORS = {
    "terrestrial": GREY,
    "semi-aquatic": SEMI,
    "non-marine aquatic": ORANGE,
    "marine": BLUE,
}
MARKERS = {
    "terrestrial": "o",
    "semi-aquatic": "^",
    "non-marine aquatic": "s",
    "marine": "o",
}
PUBLIC_PANEL = {
    "marine_model_OOF_distribution": "marine_model_OOF_distribution",
    "aquatic_v2_model_OOF_distribution": "binary_aquatic_dependence_model_OOF_distribution",
}


def ensure_dirs() -> None:
    for path in [OUTPUTS, FINAL, QUICK, QC]:
        path.mkdir(parents=True, exist_ok=True)


def load_data() -> pd.DataFrame:
    if not INPUT.exists():
        raise FileNotFoundError(INPUT)
    data = pd.read_csv(INPUT, sep="\t", keep_default_na=False)
    required = {
        "panel",
        "species_id",
        "display_name_or_species_id",
        "trait_category",
        "marine_binary",
        "aquatic_v2",
        "aquaticity_score",
        "OOF_prediction",
        "prediction_source",
        "OOF_available",
        "plotted",
        "exclusion_reason",
        "model_family",
    }
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Missing required columns: {missing}")
    return data


def write_public_data(data: pd.DataFrame) -> Path:
    out = data.copy()
    out["panel"] = out["panel"].map(PUBLIC_PANEL).fillna(out["panel"])
    out = out.rename(
        columns={
            "display_name_or_species_id": "display_name",
            "aquatic_v2": "binary_aquatic_endpoint",
        }
    )
    out["model_family"] = out["model_family"].str.replace("aquatic_v2", "binary aquatic-dependence", regex=False)
    out["exclusion_reason"] = out["exclusion_reason"].str.replace(
        "semi-aquatic aquatic_v2=0.5 excluded from accepted binomial aquatic_v2 LASSO/AUC",
        "semi-aquatic binary aquatic-dependence endpoint=0.5 excluded from accepted binomial aquatic-dependence LASSO/AUC",
        regex=False,
    )
    columns = [
        "panel",
        "species_id",
        "display_name",
        "trait_category",
        "marine_binary",
        "binary_aquatic_endpoint",
        "aquaticity_score",
        "OOF_prediction",
        "prediction_source",
        "OOF_available",
        "plotted",
        "exclusion_reason",
        "model_family",
    ]
    path = FINAL / "Figure2C_endpointfix_OOF_main_vFinal_data.tsv"
    out[columns].to_csv(path, sep="\t", index=False)
    shutil.copy2(path, OUTPUTS / path.name)
    shutil.copy2(path, QC / path.name)
    return path


def category_labels(data: pd.DataFrame, panel: str) -> list[str]:
    labels: list[str] = []
    for cat in CATEGORY_ORDER:
        rows = data[(data["panel"] == panel) & (data["trait_category"] == cat)]
        plotted = rows[rows["plotted"].astype(str).str.lower() == "true"]
        if panel == "aquatic_v2_model_OOF_distribution" and cat == "semi-aquatic":
            labels.append(f"semi-aquatic\nexcluded\n(n={len(rows)})")
        elif cat == "non-marine aquatic":
            labels.append(f"non-marine\naquatic\n(n={len(plotted)})")
        else:
            labels.append(f"{cat}\n(n={len(plotted)})")
    return labels


def plot_panel(ax: plt.Axes, data: pd.DataFrame, panel: str, title: str, ylabel: str, rng: np.random.Generator) -> None:
    ax.set_title(title, fontsize=11.8, fontweight="bold", pad=7)

    if panel == "aquatic_v2_model_OOF_distribution":
        ax.axvspan(1.58, 2.42, color="#f7f7f7", zorder=0)
        ax.text(
            2,
            0.5,
            "excluded",
            ha="center",
            va="center",
            rotation=90,
            fontsize=9.5,
            color="#777777",
            zorder=1,
        )

    for idx, cat in enumerate(CATEGORY_ORDER, start=1):
        rows = data[(data["panel"] == panel) & (data["trait_category"] == cat)]
        plotted = rows[rows["plotted"].astype(str).str.lower() == "true"].copy()
        values = pd.to_numeric(plotted["OOF_prediction"], errors="coerce").dropna().to_numpy()
        if len(values) == 0:
            continue
        ax.boxplot(
            values,
            positions=[idx],
            widths=0.52,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#222222", "linewidth": 1.15},
            boxprops={"facecolor": "white", "edgecolor": "#424242", "linewidth": 0.8},
            whiskerprops={"color": "#424242", "linewidth": 0.75},
            capprops={"color": "#424242", "linewidth": 0.75},
        )
        jitter = rng.uniform(-0.16, 0.16, len(plotted))
        ax.scatter(
            np.full(len(plotted), idx) + jitter,
            pd.to_numeric(plotted["OOF_prediction"], errors="coerce"),
            s=13,
            c=COLORS[cat],
            marker=MARKERS[cat],
            edgecolor="white",
            linewidth=0.2,
            alpha=0.72,
            zorder=3,
        )

    ax.set_xlim(0.45, 4.55)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xticks(range(1, 5))
    ax.set_xticklabels(category_labels(data, panel), fontsize=8.4)
    ax.set_ylabel(ylabel, fontsize=10.2)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=8.8)
    ax.grid(axis="y", color="#e0e0e0", linewidth=0.65)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def draw(data: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    rng = np.random.default_rng(20260522)

    # Wide, shallow panel designed for Illustrator assembly into Figure 2.
    fig = plt.figure(figsize=(9.6, 3.85), facecolor="white")
    fig.text(0.015, 0.955, "C.", fontsize=15, ha="left", va="top")
    fig.text(0.08, 0.955, "Out-of-fold prediction distributions", fontsize=14, fontweight="bold", ha="left", va="top")
    fig.text(
        0.08,
        0.897,
        "All species; predictions from genus-level leave-one-out cross-validation.",
        fontsize=9.5,
        color=TEXT_GREY,
        ha="left",
        va="top",
    )
    fig.text(
        0.565,
        0.897,
        "Semi-aquatic excluded from binary aquatic-dependence model/AUC.",
        fontsize=8.8,
        color=TEXT_GREY,
        ha="left",
        va="top",
    )

    gs = fig.add_gridspec(1, 2, left=0.085, right=0.985, bottom=0.205, top=0.79, wspace=0.27)
    ax_marine = fig.add_subplot(gs[0, 0])
    ax_aquatic = fig.add_subplot(gs[0, 1])

    plot_panel(
        ax_marine,
        data,
        "marine_model_OOF_distribution",
        "Marine model",
        "Marine-model OOF probability",
        rng,
    )
    plot_panel(
        ax_aquatic,
        data,
        "aquatic_v2_model_OOF_distribution",
        "Binary aquatic-dependence model",
        "Aquatic-dependence-model OOF probability",
        rng,
    )

    for ext in ["pdf", "png", "svg"]:
        out = FINAL / f"Figure2C_endpointfix_OOF_main_vFinal.{ext}"
        fig.savefig(out, dpi=300)
        shutil.copy2(out, OUTPUTS / out.name)
        shutil.copy2(out, QUICK / out.name)
        shutil.copy2(out, QC / out.name)
    plt.close(fig)


def write_qc_note(data_path: Path) -> None:
    note = """# Figure 2C OOF Main vFinal Visual QC

- Standalone Figure 2C panel only; Figure 2A/B were not redrawn here.
- Source values are existing all-species genus-level LOOCV out-of-fold predictions.
- No new model fitting, t-test, enrichment, or GBI analysis was run.
- Visible figure text uses OOF/out-of-fold wording and does not contain `aquatic_v2` or `full-fit`.
- Semi-aquatic taxa excluded from the binary aquatic-dependence model/AUC are retained in the data table with `prediction_source = not_available_excluded`.
- Output aspect ratio is wide and shallow for manual Illustrator assembly.

Data file:
`{data_path}`
""".format(data_path=data_path)
    path = QUICK / "Figure2C_endpointfix_OOF_main_vFinal_visual_QC.md"
    path.write_text(note, encoding="utf-8")
    shutil.copy2(path, QC / path.name)


def main() -> None:
    ensure_dirs()
    data = load_data()
    public_data = write_public_data(data)
    draw(data)
    write_qc_note(public_data)
    print(FINAL / "Figure2C_endpointfix_OOF_main_vFinal.pdf")
    print(FINAL / "Figure2C_endpointfix_OOF_main_vFinal.png")
    print(FINAL / "Figure2C_endpointfix_OOF_main_vFinal.svg")
    print(public_data)


if __name__ == "__main__":
    main()
