# Workflow overview

1. Genome/CDS alignment QC and species curation.
2. Species tree and gene-wise branch-length estimation.
3. SplitAligner branch-coordinate projection.
4. GBI matrix construction.
5. Deterministic branch-state labels and global branch-level single-gene screens.
6. Nested supervised feature-selection gLOOCV and corrected LASSO preprocessing.
7. Full-data LASSO architecture summaries.
8. Figure/source-data and supplementary-table generation.

The main validation AUCs in Fig. 2 and Fig. 5A are nested with respect to supervised
t-test/FDR feature selection and LASSO fitting, but not fully nested with respect to
global phylogenomic feature construction.
