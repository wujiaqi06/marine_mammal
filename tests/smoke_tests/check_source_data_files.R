#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value=TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[1]) else "tests/smoke_tests/check_source_data_files.R"
root <- normalizePath(file.path(dirname(this_file), "..", ".."), mustWork=TRUE)
required <- c(
  "source_data/Fig2_nested_gLOOCV/Figure2B_nested_ttest_ROC_data.tsv",
  "source_data/Fig2_nested_gLOOCV/Figure2C_nested_ttest_OOF_distribution_data.tsv",
  "source_data/Fig4_corrected_full_data_architecture/Figure4A_plot_table.tsv",
  "source_data/Fig5_sensitivity_permutation_turnover/Figure5A_nested_sensitivity_table.tsv",
  "source_data/Fig5_sensitivity_permutation_turnover/Figure5B_endpointfix_source_check.tsv",
  "source_data/Fig5_sensitivity_permutation_turnover/Figure5C_module_null_summary_2col.tsv"
)
package_root <- root
missing <- required[!file.exists(file.path(package_root, required))]
if (length(missing)) {
  writeLines(paste("Missing:", paste(missing, collapse=", ")))
  quit(status=1)
}
writeLines("PASS source data files present")
