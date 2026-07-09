#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CHECK_TARGETS=(
  "$ROOT/source_data"
  "$ROOT/README.md"
  "$ROOT/docs"
  "$ROOT/data_manifest"
)
if grep -R -n -E '17,343|1,080|98\.7%|Exact legacy permutations|intercept-only collapse|collapsed to null|near-collapse' "${CHECK_TARGETS[@]}"   --exclude='Supplementary_Table_S6_QC.tsv' --exclude='Figure5_forbidden_legacy_text_check.tsv'; then
  echo "Forbidden legacy text/value found"
  exit 1
fi
echo "PASS forbidden legacy values absent from public-facing staged assets"
