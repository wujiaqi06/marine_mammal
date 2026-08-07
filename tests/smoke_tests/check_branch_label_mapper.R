args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
script_file <- normalizePath(sub("^--file=", "", script_arg[[1L]]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
mapping_dir <- file.path(root, "metadata", "branch_label_mapping")

crosswalk <- read.delim(
  file.path(mapping_dir, "endpointfix_branch_label_crosswalk_new_to_old.tsv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = character(0)
)
stopifnot(
  nrow(crosswalk) == 601L,
  !anyDuplicated(crosswalk$new_branch_label),
  !anyDuplicated(crosswalk$old_branch_label),
  all(crosswalk$mapping_status == "PASS_EXACT_SPLIT_MATCH")
)

toy <- as.data.frame(
  matrix(c("TEST_GENE", crosswalk$new_branch_label), nrow = 1L),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
names(toy) <- c("gene", crosswalk$new_branch_label)
toy[[crosswalk$new_branch_label[[1L]]]] <- "NA"
input_file <- tempfile("branch_mapper_input_", fileext = ".tsv")
output_file <- tempfile("branch_mapper_output_", fileext = ".tsv")
write.table(toy, input_file, sep = "\t", quote = FALSE, row.names = FALSE)

mapper <- file.path(mapping_dir, "apply_branch_label_crosswalk.R")
status <- system2(
  "Rscript",
  c(mapper, paste0("--input=", input_file), paste0("--output=", output_file))
)
stopifnot(status == 0L, file.exists(output_file))

mapped <- read.delim(
  output_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  colClasses = "character",
  na.strings = character(0)
)
expected <- crosswalk[order(as.integer(sub("^B", "", crosswalk$old_branch_label))), ]
expected_values <- expected$new_branch_label
expected_values[expected$new_branch_label == crosswalk$new_branch_label[[1L]]] <- "NA"
stopifnot(
  identical(names(mapped), c("gene", paste0("B", seq_len(601L)))),
  identical(as.character(mapped[1L, paste0("B", seq_len(601L))]),
            expected_values),
  identical(mapped[[crosswalk$old_branch_label[[1L]]]][[1L]], "NA")
)

cat("PASS portable branch-label mapper: 601 exact columns; literal NA retained\n")
