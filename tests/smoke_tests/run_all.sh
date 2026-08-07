#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
Rscript "$ROOT/tests/smoke_tests/check_source_data_files.R"
Rscript "$ROOT/tests/smoke_tests/check_expected_counts.R"
bash "$ROOT/tests/smoke_tests/check_forbidden_old_values.sh"
Rscript "$ROOT/scripts/build_source_data_manifest.R"
Rscript "$ROOT/tests/smoke_tests/check_nc_source_provenance.R"
Rscript "$ROOT/tests/smoke_tests/check_branch_label_mapper.R"
Rscript "$ROOT/demo/check_nc_source_data.R" --output="$ROOT/demo/output/demo_summary.tsv"
diff -u "$ROOT/demo/expected_output/demo_summary.tsv" "$ROOT/demo/output/demo_summary.tsv"
echo "PASS all package smoke tests"
