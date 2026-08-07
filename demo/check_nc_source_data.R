#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
out_arg <- grep("^--output=", args, value = TRUE)
out_file <- if (length(out_arg)) sub("^--output=", "", out_arg[[1]]) else "demo/output/demo_summary.tsv"

args_full <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "demo/check_nc_source_data.R"
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

read_csv <- function(path) utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
fail <- function(message) {
  writeLines(paste("FAIL", message))
  quit(status = 1)
}

fig6_root <- file.path(root, "source_data", "Fig6_species_ancestral_fingerprints")
figs3_root <- file.path(root, "source_data", "FigS3_ancestor_fingerprints")

fig6a <- read_csv(file.path(fig6_root, "SourceData_Fig6A_projection_profiles.csv"))
fig6bc <- read_csv(file.path(fig6_root, "SourceData_Fig6BC_species_fingerprints_long.csv"))
s3_profiles <- read_csv(file.path(figs3_root, "SourceData_FigS3_ancestor_profile_scores.csv"))
s3_long <- read_csv(file.path(figs3_root, "SourceData_FigS3_ancestor_fingerprints_long.csv"))

portrait_species <- unique(fig6bc$species)
marine_genes <- unique(fig6bc$gene[fig6bc$model == "marine"])
aquatic_genes <- unique(fig6bc$gene[fig6bc$model == "binary_aquatic_dependence"])
ancestor_labels <- unique(s3_profiles$display_label)
s3_marine_genes <- unique(s3_long$gene[s3_long$model == "marine"])
s3_aquatic_genes <- unique(s3_long$gene[s3_long$model == "binary_aquatic_dependence"])

if (length(portrait_species) != 6L) fail("Fig. 6B/C portrait species count is not 6")
if (length(marine_genes) != 71L) fail("Fig. 6 marine predictor count is not 71")
if (length(aquatic_genes) != 101L) fail("Fig. 6 aquatic predictor count is not 101")
if (length(ancestor_labels) != 13L) fail("Fig. S3 ancestor profile count is not 13")
if (!setequal(marine_genes, s3_marine_genes)) fail("Fig. 6 and Fig. S3 marine gene sets differ")
if (!setequal(aquatic_genes, s3_aquatic_genes)) fail("Fig. 6 and Fig. S3 aquatic gene sets differ")
if (nrow(fig6a) != 22L) fail("Fig. 6A source row count is not 22")

summary <- data.frame(
  item = c(
    "Figure6A_profile_rows",
    "Figure6BC_portrait_species",
    "marine_predictors",
    "aquatic_predictors",
    "FigureS3_ancestor_profiles",
    "FigureS3_fingerprint_rows",
    "Figure6_FigureS3_gene_sets_match"
  ),
  value = c(
    nrow(fig6a),
    length(portrait_species),
    length(marine_genes),
    length(aquatic_genes),
    length(ancestor_labels),
    nrow(s3_long),
    "yes"
  ),
  stringsAsFactors = FALSE
)

dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
utils::write.table(summary, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines("PASS NC source-data checks")
