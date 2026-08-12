# Analysis update addendum

Release `v1.1-analysis-update` extends the public `v1.0-endpoint-fix` code
baseline without changing endpoint-fix scientific results. It adds the current
Fig. 6 and Supplementary Fig. S3 code, lightweight Source Data and
reviewer-facing reproducibility checks under:

- `code/12_species_ancestral_fingerprints/`
- `source_data/Fig6_species_ancestral_fingerprints/`
- `source_data/FigS3_ancestor_fingerprints/`

Fig. 6 terminal projections and gene fingerprints are descriptive corrected
full-data architecture summaries after validation. Internal branches in Fig. 6A
are terminal-only descriptive projections. Fig. S3 uses the corrected full-data
71/101 predictor architectures for internal-branch fingerprints. These layers
are intentionally kept distinct in the Source Data and documentation; their
coordinates are not silently substituted for one another.

No current scientific table was recomputed during checklist hardening. The new
demo and smoke checks read the existing reviewed Source Data and verify locked
row, species and predictor counts.
