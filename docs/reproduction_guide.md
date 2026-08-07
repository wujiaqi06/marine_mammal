# Reproduction guide

This repository contains the public code side of the marine mammal endpoint-fix
analysis package for the Nature Communications submission. Large input and
output matrices are not stored in GitHub; they are deposited in Dryad under DOI
https://doi.org/10.5061/dryad.dz08kpsd4.

## Setup

1. For the lightweight reviewer demo, run `bash demo/run_demo.sh`; the paired
   Dryad dataset is not needed.
2. For full reproduction, download and unpack the paired Dryad dataset and set
   `MARINE_MAMMAL_ENDPOINTFIX_ROOT` to the reconstructed endpoint-fix project
   root. Set `MARINE_MAMMAL_NC_OUTPUT_ROOT` to a writable output directory.
3. `config/example_paths.yaml` documents the expected path roles, but scripts
   read environment variables and do not parse the YAML file.
4. Inspect `docs/reproduce_map.tsv` to map manuscript values to source files.
5. Review `metadata/branch_label_mapping/README_branch_label_mapping.md` before
   using direct SplitAligner outputs. The public Dryad branch-coordinate and GBI
   matrices are already converted to the old downstream branch-label order used
   by the manuscript analyses.

When starting from a matrix produced directly by current SplitAligner, convert
it once before downstream analysis:

```bash
Rscript metadata/branch_label_mapping/apply_branch_label_crosswalk.R \
  --input=/path/to/current_split_aligner_matrix.tsv \
  --output=/path/to/old_label_matrix.tsv
```

Do not apply this conversion to the paired Dryad matrices; those files already
use the manuscript old-label rule and order.

Example:

```bash
export MARINE_MAMMAL_ENDPOINTFIX_ROOT=/path/to/unpacked_endpointfix_project
export MARINE_MAMMAL_NC_OUTPUT_ROOT=/path/to/reproduction_outputs
bash tests/smoke_tests/run_all.sh
```

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
- `12_NC_species_ancestral_fingerprints`: current NC Fig. 6 terminal/internal
  projections and Fig. S3 ancestral-branch fingerprint source and artwork
  scripts.

The `source_data/` directory contains lightweight real source tables used to
verify figure panels and supplementary tables. The `tests/smoke_tests/`
directory contains small checks against expected endpoint-fix values.

## Lightweight figure regeneration

These commands use bundled data only and do not rerun GBI, ASR, t-tests or
LASSO fitting:

```bash
python3 code/12_NC_species_ancestral_fingerprints/build_Figure6A_projection_profiles.py
python3 code/12_NC_species_ancestral_fingerprints/build_Figure6_layout.py
Rscript code/12_NC_species_ancestral_fingerprints/build_FigureS3_two_column.R
```

They write source-consistent draft artwork under `reproduction_outputs/`.
Final manuscript typography was adjusted in Adobe Illustrator, so raster/vector
layout details need not be byte-identical to the submitted artwork; plotted
scientific values are locked to the bundled source tables.

## Full NC Fig. 6 / Fig. S3 data regeneration

After restoring the endpoint-fix project structure from Dryad, run the scripts
below in order. They use the environment-controlled project and output roots:

```bash
Rscript code/12_NC_species_ancestral_fingerprints/run_nc_arc_terminal_only_projection.R
Rscript code/12_NC_species_ancestral_fingerprints/build_fullData_fingerprints.R
Rscript code/12_NC_species_ancestral_fingerprints/build_FigureS3_ancestor_fingerprints.R
```

The first script regenerates the terminal-only projection layer. The second
reconstructs the corrected full-data fingerprint layer after validation. The
third projects selected internal branches using the 71/101 corrected full-data
architectures and verifies coefficient reconstruction. See the code-directory
README for evidence-layer boundaries.

## Adapting to other data

SplitAligner can be run on another fixed species tree and gene trees following
`code/04_SplitAligner/README.md`. Downstream scripts require matching terminal
taxon identifiers, a gene-by-branch matrix, a trait table and branch identifiers
mapped to the endpoint-fix downstream label contract. Users must substitute
their own paths and inputs and should update expected-value tests rather than
interpreting this study's locked counts as generic software outputs.

## Provenance-script status

`provenance_scripts/nonportable_script_manifest.tsv` distinguishes older scripts
that are distributed in path-scrubbed form from historical scripts that remain
omitted. This prevents the manifest from implying that a distributed script is
absent. Portable/reviewable workflow scripts are under `code/`.
