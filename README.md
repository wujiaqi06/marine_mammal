# Marine mammal endpoint-fix reproducibility code

This repository contains custom analysis and figure-generation scripts for the
marine mammal genome evolution manuscript submitted to Nature Communications.
It contains code, lightweight source data, documentation, metadata, environment
captures, branch-label mapping metadata and smoke tests.

Large input and output data are deposited in Dryad under DOI
https://doi.org/10.5061/dryad.dz08kpsd4, including fixed-topology gene
trees/branch lengths, branch-coordinate matrices, GBI matrices, null outputs,
screening outputs and supplementary tables. The Dryad DOI is reserved; public
release is pending manuscript publication. During peer review, reviewers should
use the private Dryad link provided in the manuscript, not GitHub, to access
large data files.

SplitAligner is available separately at https://github.com/wujiaqi06/SplitAligner.
Release tag for this Nature Communications submission: `v1.0-endpoint-fix`.

## Branch-label mapping warning

SplitAligner's branch-label annotation rule changed after the downstream
manuscript analyses were built. Current direct SplitAligner outputs must be
mapped to the old downstream branch-label rule/order before reproducing the
manuscript t-test, deterministic ASR label joins, LASSO, Source Data or figure
outputs. The public Dryad branch-coordinate and GBI matrices are already
pre-converted to old labels/order. Mapping tables are provided in
`metadata/branch_label_mapping/`; the historical branch-label provenance script
is isolated under `provenance_scripts/not_portable/`.

## Minimal workflow

1. Download/unpack the Dryad dataset.
2. Copy `config/example_paths.yaml` and update the local paths.
3. Inspect `docs/reproduce_map.tsv` and `source_data/`.
4. Run smoke tests under `tests/smoke_tests/`.

See `docs/reproduction_guide.md` for a more detailed reproduction outline and
`docs/dryad_manifest.md` for large-data boundaries.

## Expected key values

- 302 terminal species.
- Fig. 2 nested AUCs: marine 0.936; binary aquatic-dependence 0.826.
- Fig. 4 corrected full-data predictors: 71 marine, 101 binary aquatic-dependence, 148 union.
- Fig. 5B positive-count-matched permutation: 894 FDR genes, 98.0% slow.

## Citation

Please cite the associated manuscript, the Dryad dataset
https://doi.org/10.5061/dryad.dz08kpsd4, and SplitAligner. The Dryad DOI is
reserved and will become publicly available when the dataset is released with
the manuscript.


## License and external data boundary

Custom code written for this study is released under the MIT license. The
external TOGA-derived mammalian coding-gene alignments from Michael Hiller Lab
resources are not redistributed here and are not relicensed by this repository.
Study-generated large processed data products are provided through the paired
Dryad dataset.

## Portable scripts and provenance scripts

Scripts intended for public reuse are kept under `code/` and should be run after
editing `config/example_paths.yaml` or setting local path environment variables.
Historical scripts that contained local absolute paths or internal package-build
state are listed in `provenance_scripts/nonportable_script_manifest.tsv` but are
not distributed in this public GitHub release. Portable/reviewable scripts are
kept under `code/`.

The `scripts/` directory is a lightweight public entry point that directs users
to the curated `code/` workflow directories.
