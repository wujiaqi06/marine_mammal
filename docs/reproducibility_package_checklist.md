# Code and software reproducibility checklist

Status refers to this review package, not to a future public tag.

| Checklist item | Status | Package evidence |
|---|---|---|
| Source code or compiled standalone software | PASS | Custom source code is under `code/`; the package is script based and needs no compilation. |
| Small real or simulated demonstration dataset | PASS | Real lightweight tables are under `source_data/`; `demo/` uses them without the large Dryad matrices. |
| Operating systems and dependencies, including versions | PASS | `README.md`, `env/R_required_packages.tsv`, `env/python_requirements.txt` and environment captures under `env/`. |
| Versions/platforms tested | PASS | macOS Sonoma 14.4 on Apple silicon; R 4.4.2, Python 3.9.6 and Perl 5.34.1. |
| Non-standard hardware | PASS | None for the demo; one CPU core and less than 1 GB RAM are sufficient. Full upstream workflows require a multicore workstation or cluster; analyses used up to 12 cores where supported. |
| Installation instructions | PASS | `README.md` and `docs/reproduction_guide.md`. |
| Typical installation time | PASS | Approximately 5-15 minutes for R/Python dependencies on a normal networked desktop; checkout/unpack is normally under 1 minute. |
| Demo run instructions | PASS | `bash demo/run_demo.sh`. |
| Demo expected output | PASS | `demo/expected_output/demo_summary.tsv` and `docs/expected_outputs.md`. |
| Typical demo run time | PASS | Less than 5 seconds after R is available on the tested desktop; exact package-preparation timing is recorded in `docs/release_command_log.tsv`. |
| Instructions for users' own data | PASS | `README.md` and the adapting-to-other-data section of `docs/reproduction_guide.md`. |
| Reproduction instructions | PASS | Stage map in `docs/reproduction_guide.md`, `docs/reproduce_map.tsv`, script READMEs and Dryad manifest. |
| License | PASS | MIT License for study-authored code. External TOGA/Hiller data are neither redistributed nor relicensed. |
| Open repository link | PASS | <https://github.com/wujiaqi06/marine_mammal>; release tag `v1.1-analysis-update`. The immutable endpoint-fix baseline remains available as `v1.0-endpoint-fix`. |
| Detailed code-function description / pseudocode location | PASS | Manuscript Methods plus `docs/workflow_overview.md`, `docs/reproduction_guide.md`, `docs/reproduce_map.tsv` and per-stage code READMEs. |

The paired Dryad DOI is <https://doi.org/10.5061/dryad.dz08kpsd4>. Public
release is pending manuscript publication; the private reviewer URL is supplied
only through the manuscript/editorial system and is intentionally excluded from
this software archive.
