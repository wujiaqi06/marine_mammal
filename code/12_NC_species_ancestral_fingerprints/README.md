# NC Fig. 6 and Supplementary Fig. S3 workflows

This directory contains the current NC-only additions to the frozen endpoint-
fix release. No scientific values were changed when the scripts were copied
into the submission package.

## Lightweight artwork regeneration

These commands use bundled source tables and write to `reproduction_outputs/`:

```bash
python3 code/12_NC_species_ancestral_fingerprints/build_Figure6A_projection_profiles.py
python3 code/12_NC_species_ancestral_fingerprints/build_Figure6_layout.py
Rscript code/12_NC_species_ancestral_fingerprints/build_FigureS3_two_column.R
```

Final typography and panel spacing were adjusted in Adobe Illustrator. The
scripts regenerate source-consistent draft artwork, not the Illustrator edit.

## Full profile reconstruction

The following scripts use the full endpoint-fix directory structure and large
matrices. Set `MARINE_MAMMAL_ENDPOINTFIX_ROOT` and optionally
`MARINE_MAMMAL_NC_OUTPUT_ROOT` before running:

```bash
Rscript code/12_NC_species_ancestral_fingerprints/run_nc_arc_terminal_only_projection.R --cores=12
Rscript code/12_NC_species_ancestral_fingerprints/build_fullData_fingerprints.R
Rscript code/12_NC_species_ancestral_fingerprints/build_FigureS3_ancestor_fingerprints.R
```

`run_nc_arc_terminal_only_projection.R` generates descriptive internal-branch
projections from models trained on terminal species. Internal branches are
projection targets only, not LASSO training or validation observations.

`build_fullData_fingerprints.R` reconstructs final full-data terminal profiles
and verifies the 71/101 selected-predictor architectures against archived
coefficients. `build_FigureS3_ancestor_fingerprints.R` applies the same final
architectures to selected internal branches and reconstructs scores including
the intercept.
