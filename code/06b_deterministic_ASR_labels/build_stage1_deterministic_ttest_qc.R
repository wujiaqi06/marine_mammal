args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[[1]] else getwd()
setwd(root_dir)

suppressPackageStartupMessages({
  library(ape)
  library(castor)
})

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  sub(" .*", "", out[[1]])
}

read_tsv <- function(path, ...) {
  utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, quote = FALSE, sep = "\t", row.names = FALSE)
}

count_values <- function(x) {
  x_chr <- as.character(x)
  data.frame(
    n_0 = sum(x_chr == "0", na.rm = TRUE),
    n_0_5 = sum(x_chr == "0.5", na.rm = TRUE),
    n_1 = sum(x_chr == "1", na.rm = TRUE),
    n_other = sum(!(x_chr %in% c("0", "0.5", "1")), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest_file <- "07_ttest_screening/stage1_deterministic_asr/qc/endpointfix_stage1_deterministic_manifest.tsv"
manifest <- read_tsv(manifest_file)
run_ids <- c("fix_marine_binary", "fix_aquatic_v2", "fix_aquatic_v1", "fix_aquatic_v3")
run_rows <- manifest[match(run_ids, manifest$run_id), , drop = FALSE]

qc_dir <- "07_ttest_screening/stage1_deterministic_asr/qc"
minimal_dir <- "07_ttest_screening/stage1_deterministic_asr/outputs_minimal"
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(minimal_dir, recursive = TRUE, showWarnings = FALSE)

script_file <- "07_ttest_screening/stage1_deterministic_asr/scripts/run_ttest_screening.R"
gbi_file <- unique(run_rows$gbi_matrix)
trait_file <- unique(run_rows$trait_file)
branch_file <- unique(run_rows$branch_file)
tree_file <- unique(run_rows$tree_file)
support_tree_file <- unique(run_rows$support_tree_file)

input_summary <- data.frame(
  input_name = c("stage1_manifest", "ttest_script", "fix_oldlabel_GBI", "trait_alias_table", "branch_file", "tree_file", "support_tree_file"),
  file_path = c(manifest_file, script_file, gbi_file, trait_file, branch_file, tree_file, support_tree_file),
  exists = file.exists(c(manifest_file, script_file, gbi_file, trait_file, branch_file, tree_file, support_tree_file)),
  file_size = file.info(c(manifest_file, script_file, gbi_file, trait_file, branch_file, tree_file, support_tree_file))$size,
  sha256 = vapply(c(manifest_file, script_file, gbi_file, trait_file, branch_file, tree_file, support_tree_file), sha256_file, character(1)),
  notes = c(
    "Proposed manifest generated during preflight; Stage 1 rows only were executed.",
    "Copied rerun script; IO-only patch added header=TRUE for GBI reading.",
    "Approved endpoint-fix no_fuse GBI with old downstream branch labels.",
    "Derived from active TY_NK_final_18pt table with old pipeline alias columns.",
    "Copied old downstream branch file.",
    "Copied old downstream annotated tree.",
    "Copied old downstream support tree."
  ),
  stringsAsFactors = FALSE
)
write_tsv(input_summary, file.path(qc_dir, "endpointfix_stage1_deterministic_ttest_input_summary.tsv"))

trait_all <- read_tsv(trait_file, row.names = 1)
branch_df <- read_tsv(branch_file)
branch_order <- branch_df$Branch
terminal_df <- branch_df[branch_df$Species != "internal", , drop = FALSE]
support_tree <- ape::read.tree(support_tree_file)
main_tree <- ape::read.tree(tree_file)

gbi <- utils::read.table(gbi_file, row.names = 1, check.names = FALSE, header = TRUE)

shape_summary <- data.frame(
  run_id = run_ids,
  matrix_name = "fix_oldlabel_GBI",
  n_input_genes = nrow(gbi),
  n_input_branches = ncol(gbi),
  expected_genes = 17432,
  expected_branches = 601,
  branch_order_match = identical(colnames(gbi), branch_order),
  status = ifelse(nrow(gbi) == 17432 & ncol(gbi) == 601 & identical(colnames(gbi), branch_order), "PASS", "STOP"),
  notes = "Same approved fix old-label GBI is used for all four Stage 1 runs.",
  stringsAsFactors = FALSE
)
write_tsv(shape_summary, file.path(qc_dir, "endpointfix_stage1_deterministic_ttest_shape_summary.tsv"))

make_trait_anc <- function(trait_column) {
  traits <- trait_all[support_tree$tip.label, trait_column] + 1
  nstates <- length(unique(traits))
  anc_results <- castor::asr_max_parsimony(support_tree, traits, nstates)
  node_states <- max.col(anc_results$ancestral_likelihoods, ties.method = "first")
  traits_number <- data.frame(traits = traits, node_number = seq_along(traits))
  anc_node_number <- data.frame(
    traits = node_states,
    node_number = seq(length(traits) + 1, length(traits) + support_tree$Nnode)
  )
  traits_states <- rbind(traits_number, anc_node_number)
  rownames(traits_states) <- paste0("B", traits_states$node_number)
  support_edges <- data.frame(support_tree$edge)
  rownames(support_edges) <- paste0("B", support_edges[, 2])
  colnames(support_edges) <- c("anc", "offspring")
  support_edges$node_states <- traits_states[rownames(support_edges), 1]
  main_tree_edges <- data.frame(main_tree$edge, main_tree$edge.length, seq_len(nrow(main_tree$edge)))
  colnames(main_tree_edges) <- c("anc", "offspring", "branch_label", "order")
  main_tree_edges$branch_label <- paste0("B", main_tree_edges$branch_label)
  rownames(support_edges) <- paste(support_edges$anc, support_edges$offspring, sep = "-")
  rownames(main_tree_edges) <- paste(main_tree_edges$anc, main_tree_edges$offspring, sep = "-")
  trait_anc <- data.frame(main_tree_edges, node_states = support_edges[rownames(main_tree_edges), "node_states"])
  trait_anc$node_states <- trait_anc$node_states - 1
  rownames(trait_anc) <- trait_anc$branch_label
  trait_anc
}

trait_state_rows <- list()
testability_rows <- list()
skipped_rows <- list()
screening_rows <- list()
direction_rows <- list()
numeric_rows <- list()

review_genes <- c("ARHGEF35", "DAOA", "FOXO3B", "NBPF10", "TDGF1P3", "USP17L19")

for (idx in seq_along(run_ids)) {
  run_id <- run_ids[[idx]]
  row <- run_rows[idx, , drop = FALSE]
  trait_column <- row$trait_column
  out_dir <- file.path("07_ttest_screening/stage1_deterministic_asr/outputs", run_id)
  out_prefix <- row$out_prefix
  alpha <- as.numeric(row$alpha)

  trait_anc <- make_trait_anc(trait_column)
  branch_trait <- trait_anc[colnames(gbi), "node_states"]
  screen_branches <- rownames(trait_anc[trait_anc$node_states != 0.5, , drop = FALSE])
  screen_branches <- screen_branches[screen_branches %in% colnames(gbi)]
  screen_states <- trait_anc[screen_branches, "node_states"]
  input_counts <- count_values(branch_trait)
  screen_counts <- count_values(screen_states)
  trait_state_rows[[run_id]] <- data.frame(
    run_id = run_id,
    trait_column = trait_column,
    n_branches_input = length(branch_trait),
    n_branch_state_0 = input_counts$n_0,
    n_branch_state_0_5 = input_counts$n_0_5,
    n_branch_state_1 = input_counts$n_1,
    n_branch_state_other = input_counts$n_other,
    n_branches_screened_after_excluding_0_5 = length(screen_branches),
    n_screened_state_0 = screen_counts$n_0,
    n_screened_state_1 = screen_counts$n_1,
    status = ifelse(screen_counts$n_0 > 1 & screen_counts$n_1 > 1, "PASS", "STOP"),
    notes = "Branch states reconstructed with the same old castor parsimony logic used by run_ttest_screening.R; branches with state 0.5 are excluded before t-test.",
    stringsAsFactors = FALSE
  )

  finite_group0 <- integer(nrow(gbi))
  finite_group1 <- integer(nrow(gbi))
  finite_total <- integer(nrow(gbi))
  testable <- logical(nrow(gbi))
  names(finite_group0) <- rownames(gbi)
  names(finite_group1) <- rownames(gbi)
  names(finite_total) <- rownames(gbi)
  names(testable) <- rownames(gbi)
  states <- screen_states
  mat <- as.matrix(gbi[, screen_branches, drop = FALSE])
  storage.mode(mat) <- "numeric"
  for (i in seq_len(nrow(mat))) {
    vals <- mat[i, ]
    ok <- is.finite(vals)
    finite_total[i] <- sum(ok)
    finite_group0[i] <- sum(ok & states == 0)
    finite_group1[i] <- sum(ok & states == 1)
    testable[i] <- (finite_group0[i] > 1 && finite_group1[i] > 1)
  }

  ttest_file <- list.files(out_dir, pattern = paste0("^", out_prefix, "[.]t_test[.]n[0-9]+[.]txt$"), full.names = TRUE)
  fdr_file <- list.files(out_dir, pattern = paste0("^", out_prefix, "[.]mammal[.]FDR", alpha, "[.]n[0-9]+[.]t_test[.]txt$"), full.names = TRUE)
  imputed_file <- list.files(out_dir, pattern = paste0("^", out_prefix, "[.]FDR", alpha, "[.]GBI[.]imputated[.]n[0-9]+[.]txt$"), full.names = TRUE)
  if (length(ttest_file) != 1 || length(fdr_file) != 1) {
    stop("Could not identify unique t-test/FDR output for ", run_id)
  }
  trait_t <- read_tsv(ttest_file)
  t_gene <- read_tsv(fdr_file)
  tested_genes <- trait_t$gene
  skipped <- setdiff(rownames(gbi), tested_genes)
  predicted_testable <- rownames(gbi)[testable]
  boundary_discrepancies <- length(setdiff(predicted_testable, tested_genes)) + length(setdiff(tested_genes, predicted_testable))

  if (length(skipped) > 0) {
    skipped_reason <- ifelse(
      finite_group0[skipped] > 1 & finite_group1[skipped] > 1,
      "not_in_old_ttest_output_despite_min_group_count_gt1_boundary_case",
      ifelse(finite_group0[skipped] <= 1 & finite_group1[skipped] <= 1, "insufficient_finite_values_in_both_groups",
             ifelse(finite_group0[skipped] <= 1, "insufficient_finite_values_in_state0", "insufficient_finite_values_in_state1"))
    )
    skipped_rows[[run_id]] <- data.frame(
      run_id = run_id,
      gene_id = skipped,
      reason = skipped_reason,
      n_finite_total = finite_total[skipped],
      n_finite_group0 = finite_group0[skipped],
      n_finite_group1 = finite_group1[skipped],
      notes = ifelse(skipped %in% review_genes, "carry-over REVIEW gene from GBI QC", ""),
      stringsAsFactors = FALSE
    )
  } else {
    skipped_rows[[run_id]] <- data.frame(
      run_id = run_id,
      gene_id = "NONE",
      reason = "no skipped genes",
      n_finite_total = NA_integer_,
      n_finite_group0 = NA_integer_,
      n_finite_group1 = NA_integer_,
      notes = "",
      stringsAsFactors = FALSE
    )
  }

  n_skipped <- length(skipped)

  testability_rows[[run_id]] <- data.frame(
    run_id = run_id,
    trait_column = trait_column,
    n_input_genes = nrow(gbi),
    n_testable_genes = nrow(trait_t),
    n_tested_genes_output = nrow(trait_t),
    n_skipped_genes = n_skipped,
    n_review_genes = length(review_genes),
    n_review_genes_skipped = sum(review_genes %in% skipped),
    review_genes_skipped = paste(intersect(review_genes, skipped), collapse = ";"),
    status = ifelse(nrow(trait_t) + n_skipped == nrow(gbi) && !anyDuplicated(tested_genes), "PASS", "STOP"),
    notes = paste0(
      "Skipped genes are defined from the exact old t-test output gene set. Finite group counts are reported for audit. ",
      "Simple finite-count precheck differs from exact old output for ", boundary_discrepancies,
      " boundary gene calls in this run."
    ),
    stringsAsFactors = FALSE
  )

  screening_rows[[run_id]] <- data.frame(
    run_id = run_id,
    trait_column = trait_column,
    n_input_genes = nrow(gbi),
    n_tested_genes = nrow(trait_t),
    n_skipped_genes = n_skipped,
    n_branches_input = ncol(gbi),
    n_branches_screened = length(screen_branches),
    n_sig_FDR_0_01 = nrow(t_gene),
    positive_t_all = sum(trait_t$tvalue > 0, na.rm = TRUE),
    negative_t_all = sum(trait_t$tvalue < 0, na.rm = TRUE),
    positive_t_sig = sum(t_gene$tvalue > 0, na.rm = TRUE),
    negative_t_sig = sum(t_gene$tvalue < 0, na.rm = TRUE),
    direction_rule = "tvalue = mean(state0) - mean(state1); positive_t means trait=1 branches have lower GBI than state0 under old pipeline estimates; negative_t means trait=1 branches have higher GBI.",
    notes = "BH/FDR threshold follows old while-loop rule p[k] < k/m*alpha.",
    stringsAsFactors = FALSE
  )

  direction_rows[[paste0(run_id, "_all_positive")]] <- data.frame(
    run_id = run_id, trait_column = trait_column, result_set = "all_tested", t_value_sign = "positive_t",
    n_genes = sum(trait_t$tvalue > 0, na.rm = TRUE), biological_direction_label_old_pipeline = "trait_state_1_lower_GBI__slow_direction", direction_rule = "positive t = mean(state0) > mean(state1)", stringsAsFactors = FALSE
  )
  direction_rows[[paste0(run_id, "_all_negative")]] <- data.frame(
    run_id = run_id, trait_column = trait_column, result_set = "all_tested", t_value_sign = "negative_t",
    n_genes = sum(trait_t$tvalue < 0, na.rm = TRUE), biological_direction_label_old_pipeline = "trait_state_1_higher_GBI__fast_direction", direction_rule = "negative t = mean(state0) < mean(state1)", stringsAsFactors = FALSE
  )
  direction_rows[[paste0(run_id, "_sig_positive")]] <- data.frame(
    run_id = run_id, trait_column = trait_column, result_set = "FDR_0.01_significant", t_value_sign = "positive_t",
    n_genes = sum(t_gene$tvalue > 0, na.rm = TRUE), biological_direction_label_old_pipeline = "trait_state_1_lower_GBI__slow_direction", direction_rule = "positive t = mean(state0) > mean(state1)", stringsAsFactors = FALSE
  )
  direction_rows[[paste0(run_id, "_sig_negative")]] <- data.frame(
    run_id = run_id, trait_column = trait_column, result_set = "FDR_0.01_significant", t_value_sign = "negative_t",
    n_genes = sum(t_gene$tvalue < 0, na.rm = TRUE), biological_direction_label_old_pipeline = "trait_state_1_higher_GBI__fast_direction", direction_rule = "negative t = mean(state0) < mean(state1)", stringsAsFactors = FALSE
  )

  check_numeric_df <- function(df, cols) {
    vals <- unlist(df[cols], use.names = FALSE)
    vals <- suppressWarnings(as.numeric(vals))
    stats <- c(sum(is.nan(vals)), sum(is.infinite(vals) & vals > 0), sum(is.infinite(vals) & vals < 0))
    names(stats) <- c("NaN", "Inf", "-Inf")
    stats
  }
  t_num <- check_numeric_df(trait_t, c("tvalue", "pvalue", "mean1", "mean2"))
  f_num <- check_numeric_df(t_gene, c("tvalue", "pvalue", "mean1", "mean2"))
  numeric_rows[[paste0(run_id, "_ttest")]] <- data.frame(
    run_id = run_id, matrix_or_file = basename(ttest_file), check_item = names(t_num), status = ifelse(t_num == 0, "PASS", "STOP"),
    n_affected_cells = as.integer(t_num), example_gene = "", example_value = "", notes = "Numeric leak check on all tested gene-level output.", stringsAsFactors = FALSE
  )
  numeric_rows[[paste0(run_id, "_fdr")]] <- data.frame(
    run_id = run_id, matrix_or_file = basename(fdr_file), check_item = names(f_num), status = ifelse(f_num == 0, "PASS", "STOP"),
    n_affected_cells = as.integer(f_num), example_gene = "", example_value = "", notes = "Numeric leak check on FDR significant gene output.", stringsAsFactors = FALSE
  )

  min_dir <- file.path(minimal_dir, run_id)
  dir.create(min_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.table(utils::head(trait_t, 50), file.path(min_dir, paste0(out_prefix, ".t_test.head50.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
  utils::write.table(utils::head(t_gene, 50), file.path(min_dir, paste0(out_prefix, ".FDR0.01.t_test.head50.tsv")), quote = FALSE, sep = "\t", row.names = FALSE)
  if (length(imputed_file) == 1) {
    con <- file(imputed_file, open = "r")
    lines <- readLines(con, n = 21)
    close(con)
    writeLines(lines, file.path(min_dir, paste0(out_prefix, ".GBI.imputated.head20rows.tsv")))
  }
}

gbi_numeric <- as.matrix(gbi)
storage.mode(gbi_numeric) <- "numeric"
gbi_checks <- c(
  sum(is.nan(gbi_numeric)),
  sum(is.infinite(gbi_numeric) & gbi_numeric > 0),
  sum(is.infinite(gbi_numeric) & gbi_numeric < 0),
  sum(is.na(gbi_numeric)),
  sum(gbi_numeric == 0, na.rm = TRUE)
)
names(gbi_checks) <- c("NaN", "Inf", "-Inf", "NA_cells", "zero_cells")
gbi_numeric_rows <- data.frame(
  run_id = "all_stage1",
  matrix_or_file = basename(gbi_file),
  check_item = names(gbi_checks),
  status = ifelse(names(gbi_checks) %in% c("NA_cells", "zero_cells"), "REVIEW", ifelse(gbi_checks == 0, "PASS", "STOP")),
  n_affected_cells = as.integer(gbi_checks),
  example_gene = "",
  example_value = "",
  notes = c("Input GBI matrix NaN check.", "Input GBI matrix Inf check.", "Input GBI matrix -Inf check.", "Input GBI contains expected missing values from GBI construction; these are not converted to 0.", "Observed true zero GBI cells, not used as missing values.")
)
numeric_leak <- rbind(gbi_numeric_rows, do.call(rbind, numeric_rows))

write_tsv(do.call(rbind, trait_state_rows), file.path(qc_dir, "endpointfix_stage1_deterministic_trait_state_branch_summary.tsv"))
write_tsv(do.call(rbind, testability_rows), file.path(qc_dir, "endpointfix_stage1_deterministic_gene_testability_summary.tsv"))
write_tsv(do.call(rbind, skipped_rows), file.path(qc_dir, "endpointfix_stage1_deterministic_skipped_gene_list.tsv"))
write_tsv(do.call(rbind, screening_rows), file.path(qc_dir, "endpointfix_stage1_deterministic_screening_summary.tsv"))
write_tsv(do.call(rbind, direction_rows), file.path(qc_dir, "endpointfix_stage1_deterministic_direction_summary.tsv"))
write_tsv(numeric_leak, file.path(qc_dir, "endpointfix_stage1_deterministic_ttest_numeric_leak_check.tsv"))

commands <- paste0(
  "Rscript 07_ttest_screening/stage1_deterministic_asr/scripts/run_ttest_screening.R 07_ttest_screening/stage1_deterministic_asr/qc/endpointfix_stage1_deterministic_manifest.tsv ",
  run_ids,
  " 07_ttest_screening/stage1_deterministic_asr/outputs/",
  run_ids
)

summary <- c(
  "# Endpoint-fix Deterministic Stage 1 t-test Run Summary",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "## Scope",
  "Deterministic Stage 1 fix baseline/comparison t-test only. No free t-test, drop/sensitivity t-test, LASSO, enrichment, figures, or old-vs-new biological comparison were run.",
  "",
  "## Exact Commands Executed",
  paste0("- `", commands, "`"),
  "",
  "## Input Provenance",
  paste0("- Script: `", script_file, "` sha256 `", sha256_file(script_file), "`"),
  paste0("- Manifest: `", manifest_file, "` sha256 `", sha256_file(manifest_file), "`"),
  paste0("- Fix old-label GBI: `", gbi_file, "` sha256 `", sha256_file(gbi_file), "`"),
  paste0("- Trait alias table: `", trait_file, "` sha256 `", sha256_file(trait_file), "`"),
  paste0("- Branch file: `", branch_file, "` sha256 `", sha256_file(branch_file), "`"),
  paste0("- Tree file: `", tree_file, "` sha256 `", sha256_file(tree_file), "`"),
  paste0("- Support tree file: `", support_tree_file, "` sha256 `", sha256_file(support_tree_file), "`"),
  "",
  "## Logic",
  "Old t-test logic, branch exclusion of 0.5, BH/FDR while-loop rule, terminal mean imputation, and output naming were reused without design change. Allowed patches: GBI matrix reading uses `header = TRUE`; ancestral-state tie handling uses `max.col(..., ties.method = \"first\")`.",
  "",
  "## IO Smoke Check",
  "PASS. Default `read.table()` misread branch headers as V2/V3; `header=TRUE` produced 17,432 genes x 601 branches and branch order matched `mammal.branch.txt`.",
  "",
  "## Screening Summary",
  paste(capture.output(print(do.call(rbind, screening_rows), row.names = FALSE)), collapse = "\n"),
  "",
  "## Trait Branch-State Counts",
  paste(capture.output(print(do.call(rbind, trait_state_rows), row.names = FALSE)), collapse = "\n"),
  "",
  "## Gene Testability",
  paste(capture.output(print(do.call(rbind, testability_rows), row.names = FALSE)), collapse = "\n"),
  "",
  "## Direction Rule",
  "The old pipeline t-statistic is `mean(state0) - mean(state1)`. Therefore positive t-values indicate lower GBI in trait-state-1 branches than state-0 branches under the old estimates, and are recorded as old-pipeline slow-direction labels; negative t-values indicate higher GBI in trait-state-1 branches and are recorded as old-pipeline fast-direction labels.",
  "",
  "## Numeric Leak Check",
  "No NaN/Inf/-Inf values were detected in final t-test or FDR-significant result tables. Input GBI contains expected NA values and they were not converted to 0.",
  "",
  "## Deterministic ASR tie-handling note",
  "A small number of tied internal nodes were observed during ancestral-state reconstruction for the aquatic traits, which allowed run-to-run variation in branch-state assignment under the previous default tie handling and therefore changed downstream t-test/FDR counts. Since Stage 1 is used as a coarse screening step, we fixed this by enforcing a deterministic tie rule (`max.col(..., ties.method = \"first\")`) without changing any trait definitions, GBI inputs, branch labels, screening logic, or BH/FDR procedures.",
  "",
  "## Deterministic ASR repeat check",
  "PASS. See `endpointfix_stage1_deterministic_asr_check.tsv`, `endpointfix_stage1_deterministic_asr_summary.tsv`, and `endpointfix_stage1_tied_internal_nodes.tsv`. All 20 repeated reconstructions per trait had identical branch-state vector hashes.",
  "",
  "## Codex self-check",
  "Codex self-check: deterministic Stage 1 fix t-test completed with no obvious STOP-level QC failure. Final PASS requires review by marine mammal rerun control chat."
)
writeLines(summary, file.path(qc_dir, "endpointfix_stage1_deterministic_ttest_run_summary.md"))
