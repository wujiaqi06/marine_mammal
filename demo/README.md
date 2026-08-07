# Reviewer demo

Run from the repository root:

```bash
bash demo/run_demo.sh
```

The demo reads only bundled real Source Data. It writes
`demo/output/demo_summary.tsv`, compares it byte-for-byte with the locked table
under `demo/expected_output/`, and then runs the endpoint-fix smoke checks.

Expected runtime is less than 5 seconds on a normal desktop after R is
installed. The demo does not download data, fit models, or require special
hardware.
