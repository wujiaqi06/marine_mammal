#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/smoke_tests/check_expected_counts.R"
pkg <- normalizePath(file.path(dirname(this_file), "..", ".."), mustWork=TRUE)
fail <- function(msg) { writeLines(msg); quit(status=1) }
trait <- file.path(pkg, "source_data/Fig1_trait_framework/trait_table.mammal302.active_TY_NK_final_18pt.tsv")
if (!file.exists(trait)) fail("Missing trait table")
t <- read.delim(trait, check.names=FALSE)
if (nrow(t) != 302) fail(paste("Expected 302 species, got", nrow(t)))
cat("PASS trait table row count = 302\n")
fig5b <- file.path(pkg, "source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_source_check.tsv")
if (!file.exists(fig5b)) fail("Missing Fig5B source check")
s <- read.delim(fig5b, check.names=FALSE)
txt <- paste(capture.output(print(s)), collapse=" ")
for (value in c("894","876","18","200")) {
  if (!grepl(value, txt, fixed=TRUE)) fail(paste("Expected Fig5B value not found:", value))
}
cat("PASS Fig5B public values found\n")
