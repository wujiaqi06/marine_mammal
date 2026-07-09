# Reproduction guide

This repository contains the public code side of the marine mammal endpoint-fix
analysis package for the Nature Communications submission. Large input and
output matrices are not stored in GitHub; they are deposited in Dryad under DOI
https://doi.org/10.5061/dryad.dz08kpsd4.

## Setup

1. Download or unpack the paired Dryad dataset.
2. Copy `config/example_paths.yaml` and edit the paths to point to the local
   Dryad data root, source-data root and large-matrix archive root.
3. Inspect `docs/reproduce_map.tsv` to map manuscript values to source files.
4. Review `metadata/branch_label_mapping/README_branch_label_mapping.md` before
   using direct SplitAligner outputs. The public Dryad branch-coordinate and GBI
   matrices are already converted to the old downstream branch-label order used
   by the manuscript analyses.

## Workflow outline

The curated code directories are under `code/`:

- `01_alignment_QC` and `01b_species_missing_gene_QC`: alignment and species
  missingness provenance.
- `02_species_tree_iqtree` and `03_paml_baseml_fixed_topology`: species-tree
  and fixed-topology branch-length provenance.
- `04_SplitAligner` and `05_GBI_construction`: branch-coordinate projection and
  GBI construction.
- `06_global_single_gene_screen`, `06b_deterministic_ASR_labels` and
  `06c_positive_count_permutation_control`: deterministic branch-state labels,
  branch-level screens and the positive-count-matched permutation control.
- `07_nested_ttest_baseline_gLOOCV` and `10_foldwise_ASR_sensitivity`: nested
  supervised feature-selection gLOOCV and fold-wise ASR sensitivity checks.
- `08_Figure4_Figure5_alignment`, `09_corrected_full_data_LASSO_architecture`
  and `11_figure_generation`: full-data architecture summaries, figure source
  tables and plotting scripts.

The `source_data/` directory contains lightweight source tables used to verify
published figure panels and supplementary tables. The `tests/smoke_tests/`
directory contains small checks against expected endpoint-fix values.

## Non-portable provenance scripts

Historical scripts with local absolute paths or package-build state are listed
in `provenance_scripts/nonportable_script_manifest.tsv` but are not distributed
in this public GitHub release. This keeps the public package portable and avoids
publishing local filesystem paths.
