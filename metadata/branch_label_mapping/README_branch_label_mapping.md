# Branch-label mapping required for direct SplitAligner output

The manuscript downstream analyses use the old downstream branch-label annotation
rule and old branch order. SplitAligner's annotation rule was updated after the
original downstream branch-label convention was established. Therefore, matrices
produced directly by the current SplitAligner annotation rule use endpoint-fix
native branch labels/order and must be mapped to the old downstream branch-label
rule/order before reproducing the manuscript t-test, deterministic ASR label
joins, LASSO, Source Data or figure outputs.

Required reproduction target: old downstream branch-label rule/order.

Readers do not need to perform this conversion for the public Dryad matrices:
the branch-coordinate and GBI matrices in Dryad are already pre-converted to old
downstream branch labels/order. This directory provides the audited mapping
tables for comparing or converting direct current-SplitAligner outputs. The
historical provenance/audit script is retained under
`provenance_scripts/not_portable/` because it records local historical paths.

To convert a direct current-SplitAligner matrix, run:

```bash
Rscript metadata/branch_label_mapping/apply_branch_label_crosswalk.R \
  --input=/path/to/current_split_aligner_matrix.tsv \
  --output=/path/to/old_label_matrix.tsv
```

The converter requires the complete, unique 601-column current-label set,
checks that every mapping row has `PASS_EXACT_SPLIT_MATCH`, writes columns in
old-label order `B1` through `B601`, and preserves literal `NA` cells. It refuses
to overwrite the input file. Dryad matrices should not be converted again.

Key files:

- `endpointfix_branch_label_crosswalk_new_to_old.tsv`: maps endpoint-fix
  SplitAligner-native branch labels to old downstream branch labels.
- `endpointfix_branch_label_crosswalk_old_to_new.tsv`: reverse mapping.
- `endpointfix_branch_label_crosswalk_summary.tsv`: QC summary showing that
  only 25 of 601 branch labels are unchanged and 576 require mapping.
- `apply_branch_label_crosswalk.R`: portable current-to-old matrix converter.
- Provenance/audit script: moved to
  `provenance_scripts/not_portable/metadata/branch_label_mapping/scripts_build_branch_label_crosswalk.py`
  because the original script records local historical paths.

Use the old-label Dryad matrices for manuscript reproduction. Use this mapping
only when starting from direct current-SplitAligner output.
