# NC software-submission QC report

Generated: 2026-08-06 17:51 JST

## Scope and status

- Baseline repository commit: `21e5e046c04e08474d21a11153262ca5778c9baa`.
- Public baseline tag: `v1.0-endpoint-fix`.
- NC release identifier: `v1.1-nc-submission`.
- Working branch used to prepare the NC release: `codex/nc-software-checklist`.
- The editor/reviewer ZIP adds current NC Fig. 6/Fig. S3 scripts and lightweight
  Source Data; these are also the additions carried by the NC release.
- No GBI, ASR, t-test/FDR, LASSO, enrichment or permutation analysis was rerun.
- Checklist verdict: **READY FOR EDITOR/REVIEWER ZIP**. This is a software-
  packaging/reproducibility verdict, not a scientific approval.

## Tested environment

- macOS Sonoma 14.4, Apple silicon (`arm64`).
- R 4.4.2.
- Python 3.9.6.
- Perl 5.34.1.
- Exact tested R/Python package versions are listed under `env/`.

## Demonstration and tests

- `bash demo/run_demo.sh`: PASS, real time 0.61 s.
- `bash tests/smoke_tests/run_all.sh`: PASS, real time 1.03 s.
- Distributed-script syntax: PASS for R parse, Python compile and Perl `-c`.
- All copied Python files compile: PASS.
- All copied/new R files parse: PASS.
- Source-data manifest: 46 files, hashes current: PASS.
- NC reviewed-source provenance: 12/12 files match locked SHA-256 values: PASS.

R emitted locale warnings because the isolated shell did not provide
`C.UTF-8`; these warnings did not change output or exit status.

## Lightweight figure-script execution

Scripts were invoked from outside the repository to test location-independent
path resolution:

- Fig. 6A paired profile draft: PASS, 1.15 s.
- Fig. 6 combined draft: PASS, 7.37 s.
- Fig. S3 artwork variants: PASS, 6.31 s.

The restricted verification shell required a writable Matplotlib font cache;
the documented `MPLCONFIGDIR` workaround was used. A normal user-writable
Matplotlib cache requires no workaround.

## Scientific-source preservation

The 12 NC source files under the Fig. 6 and Fig. S3 Source Data directories are
byte-identical to reviewed `analysis_nc_arc` outputs. Their relative source
paths and hashes are recorded in
`data_manifest/NC_source_data_provenance.tsv`. Renaming only reflects the
current manuscript figure numbers.

Fig. 6A internal rows and Fig. S3 internal-branch fingerprints remain distinct
evidence layers. The package documents rather than erases this distinction.

## Safety and portability checks

- Local absolute path grep over distributable code/data/docs: PASS, zero hits.
- Forbidden legacy-value/claim grep: PASS, zero hits in public-facing targets.
- Files larger than 100 MB outside `.git`: zero.
- Private Dryad reviewer URL: absent by design.
- `git diff --check`: PASS.
- Branch-label crosswalk, warning and converted-matrix contract: retained; the
  portable converter passed a complete 601-column mapping test with literal
  `NA` preservation.
- Provenance manifest inconsistency repaired: distributed scripts are no longer
  described as absent.

## Remaining author-controlled release actions

1. Release the reserved Dryad record at the manuscript-publication stage.

The private Dryad reviewer link remains in the manuscript or editorial system
only and is not stored in the public repository.

## GitHub v1.1 release-candidate refresh

Rechecked on 2026-08-07 JST before the author-controlled remote release:

- `bash demo/run_demo.sh`: PASS, real time 0.65 s.
- `bash tests/smoke_tests/run_all.sh`: PASS, real time 1.12 s.
- Source-data manifest: 46/46 files current.
- NC reviewed-source provenance: 12/12 SHA-256 values matched.
- Portable branch-label converter: all 601 columns mapped to the manuscript
  old-label order; literal `NA` values retained.
- Distributed syntax: zero R, Python or Perl failures.
- Local absolute paths, credential patterns and private reviewer URLs: zero.
- Files larger than 100 MB outside `.git`: zero.
- `git diff --check`: PASS.
- The Dryad public record remains embargoed/private during peer review; no
  private reviewer link is included in the repository.
