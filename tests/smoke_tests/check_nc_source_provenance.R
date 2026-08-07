#!/usr/bin/env Rscript

args <- commandArgs(FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "tests/smoke_tests/check_nc_source_provenance.R"
root <- normalizePath(file.path(dirname(this_file), "..", ".."), mustWork = TRUE)

if (!requireNamespace("digest", quietly = TRUE)) stop("digest is required")
manifest <- utils::read.delim(
  file.path(root, "data_manifest", "NC_source_data_provenance.tsv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

paths <- file.path(root, manifest$package_path)
if (any(!file.exists(paths))) {
  stop("Missing NC source files: ", paste(manifest$package_path[!file.exists(paths)], collapse = ", "))
}
observed <- vapply(paths, digest::digest, character(1), algo = "sha256", file = TRUE)
if (!identical(unname(observed), unname(manifest$sha256))) {
  stop("NC source-file SHA-256 mismatch: ", paste(manifest$package_path[observed != manifest$sha256], collapse = ", "))
}
if (!all(manifest$byte_identity == "PASS")) stop("NC provenance ledger contains a non-PASS byte-identity row")
cat("PASS 12 NC scientific source files match locked reviewed hashes\n")
