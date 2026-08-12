#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT_DIR="${1:-$ROOT/demo/output}"
mkdir -p "$OUT_DIR"

Rscript "$ROOT/demo/check_extended_source_data.R" --output="$OUT_DIR/demo_summary.tsv"
diff -u "$ROOT/demo/expected_output/demo_summary.tsv" "$OUT_DIR/demo_summary.tsv"

Rscript "$ROOT/tests/smoke_tests/check_source_data_files.R"
Rscript "$ROOT/tests/smoke_tests/check_expected_counts.R"
bash "$ROOT/tests/smoke_tests/check_forbidden_old_values.sh"

echo "PASS demo output matches locked expected output"
