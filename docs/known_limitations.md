# Known limitations and release boundaries

- Large study-generated matrices are deposited in Dryad rather than GitHub.
- External TOGA/Hiller raw alignment data are not redistributed or relicensed by this package.
- Direct current-SplitAligner branch labels must be mapped to the old downstream branch-label rule/order before reproducing the manuscript analyses; the Dryad matrices are already converted.
- Historical-script disposition is recorded in `provenance_scripts/nonportable_script_manifest.tsv`; the ledger distinguishes omitted scripts from path-scrubbed or provenance-labelled scripts that are distributed under `code/`.
- Dryad DOI is reserved as https://doi.org/10.5061/dryad.dz08kpsd4; Dryad public
  release is pending manuscript publication, and reviewers should use the private
  Dryad link provided in the manuscript during peer review.
- GitHub tag `v1.1-nc-submission` contains the NC Fig. 6/Fig. S3 code, lightweight Source Data and reproducibility checks; `v1.0-endpoint-fix` remains the immutable baseline.
- Lightweight figure scripts regenerate source-consistent drafts. Final manuscript typography and panel spacing were adjusted in Adobe Illustrator.
- Internal-branch projections are descriptive genomic-profile projections and not ancestral habitat assignments or held-out validation estimates.
