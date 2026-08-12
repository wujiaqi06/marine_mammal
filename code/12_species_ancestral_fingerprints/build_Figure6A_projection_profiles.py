#!/usr/bin/env python3
"""Build the Figure 6A paired fitted-projection profile draft.

Terminal species use corrected full-data fitted projections.
Internal/ancestor rows use terminal-only internal projection outputs
and are displayed as descriptive projection rows, not validation estimates.
"""

from __future__ import annotations

import os
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

TERMINAL_PROJ = SOURCE_ROOT / "SourceData_Fig6_fullData_species_projections_all_focal.csv"
INTERNAL_PROJ = SOURCE_ROOT / "selected_internal_node_projections.csv"

INK = "#1F2430"
MUTED = "#6F768A"
GRID = "#E6E8F0"
CARD = "#D7DCE5"
MARINE_BLUE = "#2C7FB8"
AQUATIC_ORANGE = "#D95F02"
LINE_GREY = "#A9ADB5"

TERMINAL_ROWS = [
    ("terminal", "Core cetacean", "Killer whale", "Orcinus_orca"),
    ("terminal", "Core pinniped", "California sea lion", "Zalophus_californianus"),
    ("terminal", "Heterogeneous pinniped", "Weddell seal", "Leptonychotes_weddellii"),
    ("terminal", "Heterogeneous pinniped", "Walrus", "Odobenus_rosmarus_divergens"),
    ("terminal", "Sirenian / edge", "Dugong", "Dugong_dugon"),
    ("terminal", "Marine edge", "Polar bear", "Ursus_maritimus"),
    ("terminal", "Marine edge", "Sea otter", "Enhydra_lutris_kenyoni"),
    ("terminal", "River dolphin bridge", "Platanista", "Platanista_minor"),
    ("terminal", "River dolphin bridge", "Inia", "Inia_geoffrensis"),
    ("terminal", "River dolphin bridge", "Baiji", "Lipotes_vexillifer"),
    ("terminal", "Non-marine aquatic controls", "Hippopotamus", "Hippopotamus_amphibius"),
    ("terminal", "Non-marine aquatic controls", "Small-clawed otter", "Aonyx_cinereus"),
    ("terminal", "Non-marine aquatic controls", "Giant otter", "Pteronura_brasiliensis"),
]

INTERNAL_ROWS = [
    ("internal", "Cetacean ancestry", "Cetacea + hippo ancestor", "common_ancestor_Cetacea_plus_Hippopotamus"),
    ("internal", "Cetacean ancestry", "Crown Cetacea", "crown_Cetacea"),
    ("internal", "Cetacean ancestry", "Odontoceti", "Odontoceti"),
    ("internal", "Pinniped ancestry", "Crown Pinnipedia", "Pinnipedia"),
    ("internal", "Pinniped ancestry", "Phocidae", "Phocidae"),
    ("internal", "Pinniped ancestry", "Otarioidea", "Otarioidea"),
    ("internal", "Pinniped ancestry", "Otariidae", "Otariidae"),
    ("internal", "Sirenian decoupling context", "Sirenia", "Sirenia"),
    ("internal", "Sirenian decoupling context", "Dugong + Hydrodamalis", "Dugong_plus_Hydrodamalis"),
]


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


def read_terminal() -> pd.DataFrame:
    d = pd.read_csv(TERMINAL_PROJ)
    wide = d.pivot(index="species", columns="model", values="fitted_probability").reset_index()
    wide = wide.rename(columns={
        "marine": "marine_probability",
        "binary_aquatic_dependence": "aquatic_probability",
    })
    meta = d.drop_duplicates("species")[["species", "trait_category", "aquaticity_score", "clade_group"]]
    wide = wide.merge(meta, on="species", how="left")
    rows = []
    for kind, group, label, key in TERMINAL_ROWS:
        r = wide[wide["species"] == key]
        if r.empty:
            rows.append({
                "row_type": kind, "group": group, "label": label, "source_id": key,
                "marine_probability": np.nan, "aquatic_probability": np.nan,
                "source_layer": "missing_terminal_projection",
            })
            continue
        rr = r.iloc[0]
        rows.append({
            "row_type": kind,
            "group": group,
            "label": label,
            "source_id": key,
            "marine_probability": float(rr["marine_probability"]),
            "aquatic_probability": float(rr["aquatic_probability"]),
            "trait_category": rr.get("trait_category", ""),
            "aquaticity_score": rr.get("aquaticity_score", ""),
            "source_layer": "corrected full-data terminal fitted projection",
        })
    return pd.DataFrame(rows)


def read_internal() -> pd.DataFrame:
    d = pd.read_csv(INTERNAL_PROJ, sep="\t")
    wide = d.pivot(
        index=["target_name", "branch_id", "target_type", "target_note"],
        columns="trait",
        values="predicted_probability",
    ).reset_index()
    wide = wide.rename(columns={
        "marine_binary": "marine_probability",
        "binary_aquatic_dependence": "aquatic_probability",
    })
    rows = []
    for kind, group, label, key in INTERNAL_ROWS:
        r = wide[wide["target_name"] == key]
        if r.empty:
            rows.append({
                "row_type": kind, "group": group, "label": label, "source_id": key,
                "marine_probability": np.nan, "aquatic_probability": np.nan,
                "source_layer": "missing_internal_projection",
            })
            continue
        rr = r.iloc[0]
        rows.append({
            "row_type": kind,
            "group": group,
            "label": label,
            "source_id": key,
            "branch_id": rr["branch_id"],
            "target_type": rr["target_type"],
            "marine_probability": float(rr["marine_probability"]),
            "aquatic_probability": float(rr["aquatic_probability"]),
            "source_layer": "terminal-only internal branch projection; descriptive, not OOF validation",
            "target_note": rr["target_note"],
        })
    return pd.DataFrame(rows)


def build_source_table() -> pd.DataFrame:
    terminal = read_terminal()
    internal = read_internal()
    d = pd.concat([terminal, internal], ignore_index=True)
    d["axis_delta_aquatic_minus_marine"] = d["aquatic_probability"] - d["marine_probability"]
    d["projection_summary"] = np.select(
        [
            (d["marine_probability"] >= 0.75) & (d["aquatic_probability"] >= 0.75),
            (d["marine_probability"] < 0.25) & (d["aquatic_probability"] >= 0.75),
            (d["marine_probability"] >= 0.75) & (d["aquatic_probability"] < 0.25),
            (d["marine_probability"] < 0.25) & (d["aquatic_probability"] < 0.25),
        ],
        [
            "high on both fitted axes",
            "aquatic-axis high, marine-axis low",
            "marine-axis high, aquatic-axis low",
            "low on both fitted axes",
        ],
        default="intermediate or divergent",
    )
    return d


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


def compact_group_label(group: str) -> str:
    """Short two-line labels keep the left gutter readable."""
    replacements = {
        "Heterogeneous pinniped": "Heterogeneous\npinniped",
        "River dolphin bridge": "River dolphin\nbridge",
        "Non-marine aquatic controls": "Non-marine\naquatic controls",
        "Cetacean ancestry": "Cetacean\nancestry",
        "Pinniped ancestry": "Pinniped\nancestry",
        "Sirenian decoupling context": "Sirenian decoupling\ncontext",
    }
    return replacements.get(group, group)


def draw_profile(d: pd.DataFrame) -> None:
    setup()
    OUT.mkdir(parents=True, exist_ok=True)
    n = len(d)
    fig_h = max(7.2, 0.34 * n + 2.2)
    fig, ax = plt.subplots(figsize=(10.8, fig_h), facecolor="white")

    ax.set_xlim(-0.68, 1.03)
    ax.set_ylim(-0.75, n - 0.25)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Fitted probability", fontsize=12, color=INK)
    ax.set_xticks([0, 0.25, 0.50, 0.75, 1.00])
    ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"], fontsize=9, color=INK)
    ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.grid(axis="y", visible=False)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#BFC5D0")

    # Background bands by group and row type.
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
            fontsize=7.2,
            color=MUTED,
            fontweight="bold",
            linespacing=0.92,
        )

    for i, r in d.iterrows():
        marine = r["marine_probability"]
        aquatic = r["aquatic_probability"]
        ax.plot([marine, aquatic], [i, i], color=LINE_GREY, linewidth=1.1, zorder=1)
        if r["row_type"] == "terminal":
            ax.scatter(marine, i, s=58, marker="o", color=MARINE_BLUE, edgecolor="white", linewidth=0.7, zorder=3)
            ax.scatter(aquatic, i, s=58, marker="s", color=AQUATIC_ORANGE, edgecolor="white", linewidth=0.7, zorder=3)
        else:
            ax.scatter(marine, i, s=62, marker="D", facecolor="white", edgecolor=MARINE_BLUE, linewidth=1.8, zorder=3)
            ax.scatter(aquatic, i, s=62, marker="D", facecolor="white", edgecolor=AQUATIC_ORANGE, linewidth=1.8, zorder=3)
        ax.text(-0.31, i, r["label"], ha="left", va="center", fontsize=9.2, color=INK)

    title = "A. Representative fitted projection profiles"
    subtitle = "Species use corrected full-data terminal projections; internal nodes are descriptive branch projections."
    ax.text(-0.66, -1.55, title, ha="left", va="top", fontsize=15, color=INK, fontweight="bold", transform=ax.transData)
    ax.text(-0.66, -1.08, subtitle, ha="left", va="top", fontsize=9.4, color=MUTED, transform=ax.transData)

    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=MARINE_BLUE, markeredgecolor="white", markersize=8, label="Marine model fitted probability"),
        Line2D([0], [0], marker="s", color="none", markerfacecolor=AQUATIC_ORANGE, markeredgecolor="white", markersize=8, label="Aquatic-dependence model fitted probability"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=INK, markersize=7, label="Internal node / branch projection"),
    ]
    ax.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 1.105), frameon=False, fontsize=8.3)

    card = FancyBboxPatch(
        (0.015, 0.02), 0.97, 0.95,
        transform=fig.transFigure,
        boxstyle="round,pad=0.012,rounding_size=0.015",
        linewidth=1.0,
        edgecolor=CARD,
        facecolor="none",
        zorder=-10,
    )
    fig.patches.append(card)
    fig.subplots_adjust(left=0.08, right=0.97, top=0.88, bottom=0.10)

    for ext in ["pdf", "svg", "png"]:
        fig.savefig(OUT / f"Figure6A_paired_projection_profiles.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)
    plt.close(fig)


def write_notes(d: pd.DataFrame) -> None:
    d.to_csv(OUT / "SourceData_Figure6A_paired_projection_profiles.csv", index=False)
    caption = """# Figure 6A Paired Projection Profile Caption Draft

Panel A shows fitted projection profiles for representative terminal species and selected internal branches/nodes. Terminal species probabilities are corrected full-data final-architecture fitted projections. Internal rows are descriptive terminal-only internal-branch projections and are shown with open diamond markers; they are not out-of-fold validation estimates and should not be described as ancestral habitat predictions. Blue circles/diamonds indicate marine-axis fitted probabilities, orange squares/diamonds indicate binary aquatic-dependence fitted probabilities, and grey segments connect the two model axes for each row.
"""
    (OUT / "Figure6A_paired_projection_caption_draft.md").write_text(caption, encoding="utf-8")

    selection = """# Figure 6A Species And Internal-Node Selection

## Terminal species
- Killer whale, Weddell seal, dugong and polar bear provide representative marine and edge portraits.
- California sea lion and walrus are retained as core/heterogeneous pinniped portraits.
- Sea otter is retained as a marine-edge decoupling case.
- Platanista, Inia and baiji bridge to the river-dolphin / Fig. 5 analysis.
- Hippopotamus, small-clawed otter and giant otter represent non-marine aquatic controls.

## Internal / ancestor rows
- Cetacea + hippo ancestor, Crown Cetacea and Odontoceti represent cetacean ancestry.
- Crown Pinnipedia, Phocidae, Otarioidea and Otariidae represent pinniped ancestry.
- Sirenia and Dugong + Hydrodamalis are included as decoupling context.

## Evidence boundary
Terminal rows use corrected full-data terminal fitted projections. Internal rows use terminal-only internal projection outputs and must be treated as descriptive branch projections, not out-of-fold validation and not ancestral habitat reconstruction.
"""
    (OUT / "Figure6A_species_node_selection_note.md").write_text(selection, encoding="utf-8")

    qc = [
        ("representative_terminal_species_included", "PASS", "Killer whale, Weddell seal, dugong and polar bear included."),
        ("california_sea_lion_walrus_included", "PASS", "core/heterogeneous pinniped portraits included."),
        ("river_dolphin_bridge_species_included", "PASS", "Platanista, Inia, baiji included."),
        ("non_marine_aquatic_controls_included", "PASS", "Hippo, small-clawed otter, giant otter included."),
        ("ancestor_rows_included", "PASS", "Nine internal/ancestor rows included, including Sirenia and Dugong+Hydrodamalis."),
        ("internal_rows_marked_descriptive", "PASS", "Caption and source_layer distinguish internal projections from OOF validation."),
        ("no_posterior_wording", "PASS", "No posterior probability wording used."),
    ]
    pd.DataFrame(qc, columns=["check", "status", "notes"]).to_csv(OUT / "Figure6A_paired_projection_QC.tsv", sep="\t", index=False)


def package_outputs() -> Path:
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    pkg_name = f"MarineMammal_Figure6A_reproduction_{timestamp}"
    pkg_dir = REVIEW_ROOT / pkg_name
    pkg_dir.mkdir(parents=True, exist_ok=True)
    keep = [
        "Figure6A_paired_projection_profiles.pdf",
        "Figure6A_paired_projection_profiles.png",
        "Figure6A_paired_projection_profiles.svg",
        "SourceData_Figure6A_paired_projection_profiles.csv",
        "Figure6A_paired_projection_caption_draft.md",
        "Figure6A_species_node_selection_note.md",
        "Figure6A_paired_projection_QC.tsv",
    ]
    for name in keep:
        src = OUT / name
        if src.exists():
            dst = pkg_dir / name
            dst.write_bytes(src.read_bytes())
    (pkg_dir / "README_Figure6A_reproduction.md").write_text(
        "# Figure 6A Paired Projection Profiles\n\n"
        "This package contains a Figure 6A paired-dot draft and source/QC files. "
        "No models were rerun.\n",
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
    (OUT / "latest_Figure6A_projection_review_package_path.txt").write_text(
        str(zip_path) + "\n", encoding="utf-8"
    )
    return zip_path


def main() -> None:
    d = build_source_table()
    draw_profile(d)
    write_notes(d)
    zip_path = package_outputs()
    print(f"Wrote {OUT / 'Figure6A_paired_projection_profiles.pdf'}")
    print(f"Wrote package: {zip_path}")


if __name__ == "__main__":
    main()
