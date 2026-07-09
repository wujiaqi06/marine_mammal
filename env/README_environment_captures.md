# Environment captures

`R_sessionInfo_current_local.txt` and `R_package_version_summary.tsv` record the current local R environment available during public-release packaging QC. They are included to make the public package more transparent, but they are not a claim that all upstream genome/tree/BaseML analyses were run in this exact local session. Upstream command-line tools and provenance are documented separately in Methods source-lock files and workflow notes.

`python_pip_freeze.txt` records the Python environment captured during package preparation after local file-url entries were removed.
