# Species and ancestral fingerprint code porting changes

The Fig. 6/Fig. S3 scripts were copied from the reviewed `reproduction_outputs`
work area and changed only for public-package portability and current output
naming. The changes were:

- replace local absolute project/output paths with documented environment
  variables;
- derive lightweight Source Data paths from the repository root;
- correct the Fig. S3 artwork script's repository-root traversal;
- avoid silent duplication of fitted-probability columns and assert agreement
  between the bundled fingerprint and projection tables;
- use current Fig. 6/Fig. S3 output names and remove internal task-number wording;
- direct figure output to ignored `reproduction_outputs/` directories.

Two pre-existing historical helper-script syntax defects were also repaired so
that the distributed source tree passes language syntax checks: an orphan final
`)` was removed from `code/01b_species_missing_gene_QC/Data.R`, and the current
sequence array name in `code/01_alignment_QC/RandomSamplingSites.pl` was
corrected from the undeclared `@cds_seq` to the already populated `@sites`.
These are mechanical source repairs; neither helper was executed to regenerate
scientific data for this package.

No coefficients, selected predictor sets, GBI values, branch IDs, profile
coordinates or plotted Source Data were edited. Byte identity of the 12 copied
reviewed scientific source files is recorded in
`data_manifest/species_ancestral_source_data_provenance.tsv`.
