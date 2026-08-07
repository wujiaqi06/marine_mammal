# Expected output locks

## Reviewer demo

Run:

```bash
bash demo/run_demo.sh
```

Expected final lines:

```text
PASS NC source-data checks
PASS source data files present
PASS trait table row count = 302
PASS Fig5B public values found
PASS forbidden legacy values absent from public-facing staged assets
PASS demo output matches locked expected output
```

The generated `demo/output/demo_summary.tsv` must be byte-identical to
`demo/expected_output/demo_summary.tsv`. The expected NC-specific locks are:

- Fig. 6A source rows: 22.
- Fig. 6B/C terminal portraits: 6.
- Marine predictors per fingerprint: 71.
- Binary aquatic-dependence predictors per fingerprint: 101.
- Fig. S3 ancestral/internal profiles: 13.
- Fig. S3 gene-contribution rows: 2,236.
- Fig. 6 and Fig. S3 gene sets: identical for each model.

See `docs/reproduce_map.tsv` for endpoint-fix manuscript values and source
files. Smoke tests also check species counts, screening counts, nested AUCs,
full-data predictor counts, Fig. 5B permutation values and forbidden legacy
values.
