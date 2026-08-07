#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(prefix) {
  hit <- args[startsWith(args, prefix)]
  if (length(hit) != 1L) {
    stop("Expected exactly one ", prefix, " argument.", call. = FALSE)
  }
  sub(prefix, "", hit, fixed = TRUE)
}

input_file <- normalizePath(get_arg("--input="), mustWork = TRUE)
output_file <- get_arg("--output=")
if (file.exists(output_file)) {
  output_file <- normalizePath(output_file, mustWork = TRUE)
}
if (identical(input_file, output_file)) {
  stop("Input and output must be different files.", call. = FALSE)
}

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
mapping_dir <- dirname(script_file)
crosswalk_file <- file.path(
  mapping_dir,
  "endpointfix_branch_label_crosswalk_new_to_old.tsv"
)

crosswalk <- read.delim(
  crosswalk_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character(0)
)
required_fields <- c("new_branch_label", "old_branch_label", "mapping_status")
if (!all(required_fields %in% names(crosswalk))) {
  stop("Crosswalk is missing required fields.", call. = FALSE)
}
if (nrow(crosswalk) != 601L ||
    anyDuplicated(crosswalk$new_branch_label) ||
    anyDuplicated(crosswalk$old_branch_label) ||
    !all(crosswalk$mapping_status == "PASS_EXACT_SPLIT_MATCH")) {
  stop("Crosswalk failed the locked 601-branch one-to-one QC contract.", call. = FALSE)
}

matrix_in <- read.delim(
  input_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = "character",
  na.strings = character(0)
)
missing_branches <- setdiff(crosswalk$new_branch_label, names(matrix_in))
if (length(missing_branches) > 0L) {
  stop(
    "Input is missing current-SplitAligner branch columns: ",
    paste(head(missing_branches, 10L), collapse = ", "),
    if (length(missing_branches) > 10L) " ..." else "",
    call. = FALSE
  )
}

id_columns <- setdiff(names(matrix_in), crosswalk$new_branch_label)
old_number <- as.integer(sub("^B", "", crosswalk$old_branch_label))
if (anyNA(old_number) || !setequal(old_number, seq_len(601L))) {
  stop("Old branch labels are not the complete B1-B601 set.", call. = FALSE)
}
crosswalk <- crosswalk[order(old_number), , drop = FALSE]

matrix_out <- matrix_in[c(id_columns, crosswalk$new_branch_label)]
names(matrix_out) <- c(id_columns, crosswalk$old_branch_label)
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
write.table(
  matrix_out,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = "NA"
)

message(
  "PASS mapped ", nrow(matrix_out), " rows and 601 branch columns to ",
  "the manuscript old-label order: ", output_file
)
