# Dryad manifest and large-data boundary

Large files required to reproduce the complete marine mammal endpoint-fix
analysis are deposited in Dryad under DOI
https://doi.org/10.5061/dryad.dz08kpsd4. Public release of the Dryad dataset is
pending manuscript publication; during peer review, use the private Dryad link
provided in the manuscript.

GitHub intentionally excludes:

- raw external TOGA/Hiller mammalian coding-gene alignment bundles;
- full fixed-topology BaseML/PAML working directories;
- full branch-coordinate matrices;
- full GBI matrices;
- compressed large matrix archives and other files larger than 100 MB.

The large-data manifest is mirrored in `data_manifest/Dryad_large_data_manifest.tsv`.
Lightweight source tables for manuscript figures and supplementary tables are
included in `source_data/`.

Branch-label mapping metadata is included under `metadata/branch_label_mapping/`.
The Dryad branch-coordinate and GBI matrices are already converted to the old
downstream branch-label order used by the manuscript analyses; readers do not
need to run branch-label conversion unless starting from direct current
SplitAligner outputs.
