#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
write_mode <- "--write" %in% args

args_full <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/build_source_data_manifest.R"
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest R package is required to verify SHA-256 values.")
}

files <- list.files(file.path(root, "source_data"), recursive = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir == FALSE]
files <- files[basename(files) != "README.md"]
files <- sort(files)

rel <- sub(paste0("^", root, "/"), "", files, fixed = FALSE)
is_extended <- grepl("^source_data/(Fig6_species_ancestral_fingerprints|FigS3_ancestor_fingerprints)/", rel)
notes <- ifelse(
  is_extended,
  "Current extended lightweight Source Data copied from reviewed Fig. 6/Fig. S3 analysis outputs; scientific values unchanged",
  "Lightweight source data mirrored from Dryad package; local source paths scrubbed to public tokens where present"
)

manifest <- data.frame(
  file_path = rel,
  size_bytes = file.info(files)$size,
  sha256 = vapply(files, digest::digest, character(1), algo = "sha256", file = TRUE),
  notes = notes,
  stringsAsFactors = FALSE
)

manifest_file <- file.path(root, "data_manifest", "source_data_manifest.tsv")
if (write_mode) {
  utils::write.table(manifest, manifest_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Wrote", nrow(manifest), "source-data manifest rows\n")
} else {
  if (!file.exists(manifest_file)) stop("Missing source-data manifest")
  locked <- utils::read.delim(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  as_text <- function(x) as.data.frame(lapply(x, as.character), stringsAsFactors = FALSE)
  if (!identical(as_text(locked), as_text(manifest))) {
    stop("Source-data manifest differs from current files; rerun with --write after reviewing changes.")
  }
  cat("PASS source-data manifest matches", nrow(manifest), "files\n")
}
