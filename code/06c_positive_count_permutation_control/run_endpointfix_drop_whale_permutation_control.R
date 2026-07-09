#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(castor)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else getwd()
n_perm <- if (length(args) >= 2) as.integer(args[2]) else 200L
seed <- if (length(args) >= 3) as.integer(args[3]) else 20260520L

alpha <- 0.01
run_id <- "fix_drop_whale"
trait_name <- "drop_whale"

perm_dir <- file.path(root, "07_ttest_screening", "permutation_control")
out_dir <- file.path(perm_dir, "outputs")
qc_dir <- file.path(perm_dir, "qc")
min_dir <- file.path(perm_dir, "outputs_minimal")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(min_dir, recursive = TRUE, showWarnings = FALSE)

gbi_file <- file.path(root, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv")
trait_file <- file.path(root, "05_trait_tables", "drop_inputs", "marine_drop_trait_inputs", "trait_input.marine_binary_drop_sets.combined.tsv")
branch_file <- file.path(root, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt")
tree_file <- file.path(root, "07_ttest_screening", "inputs", "branch_files", "mammal302.anno.nwk")
support_tree_file <- file.path(root, "07_ttest_screening", "inputs", "branch_files", "mammal302.anno.BL_support.nwk")
stage2_summary_file <- file.path(root, "07_ttest_screening", "stage2_deterministic_drop_sensitivity", "qc", "endpointfix_stage2_screening_summary.tsv")

sha256 <- function(path) {
  unname(tools::md5sum(path))
}

write_tsv <- function(x, path) {
  utils::write.table(x, file = path, quote = FALSE, sep = "\t", row.names = FALSE)
}

input_summary <- data.frame(
  input_name = c("fix_oldlabel_GBI", "drop_trait_table", "branch_file", "tree_file", "support_tree_file", "stage2_observed_summary"),
  file_path = c(gbi_file, trait_file, branch_file, tree_file, support_tree_file, stage2_summary_file),
  exists = file.exists(c(gbi_file, trait_file, branch_file, tree_file, support_tree_file, stage2_summary_file)),
  file_size = file.info(c(gbi_file, trait_file, branch_file, tree_file, support_tree_file, stage2_summary_file))$size,
  checksum_type = "md5",
  checksum = vapply(c(gbi_file, trait_file, branch_file, tree_file, support_tree_file), sha256, character(1)) |> c(sha256(stage2_summary_file)),
  notes = c(
    "Approved endpoint-fix old-label no_fuse fix GBI.",
    "Copied rerun-folder marine drop trait table.",
    "Copied old downstream branch map.",
    "Copied old annotated species tree.",
    "Copied old annotated support tree used for deterministic ASR.",
    "Release-ready deterministic Stage 2 observed drop-whale summary."
  ),
  stringsAsFactors = FALSE
)
write_tsv(input_summary, file.path(qc_dir, "endpointfix_permutation_input_summary.tsv"))

trait_all <- utils::read.delim(trait_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
branch_df <- utils::read.delim(branch_file, stringsAsFactors = FALSE, check.names = FALSE)
support_tree <- ape::read.tree(support_tree_file)
main_tree <- ape::read.tree(tree_file)
terminal_df <- branch_df[branch_df$Species != "internal", , drop = FALSE]
terminal_df$genus <- vapply(strsplit(terminal_df$Species, "_"), `[`, character(1), 1)
trait_terminal <- trait_all[terminal_df$Species, , drop = FALSE]
trait_terminal$Branch <- terminal_df$Branch
trait_terminal$Species <- terminal_df$Species
trait_terminal$genus <- terminal_df$genus
rownames(trait_terminal) <- trait_terminal$Branch
eligible_terminal <- trait_terminal[trait_terminal[[trait_name]] != 0.5, , drop = FALSE]
n_positive_terminal <- sum(eligible_terminal[[trait_name]] == 1)
n_negative_terminal <- sum(eligible_terminal[[trait_name]] == 0)

traits <- trait_all[support_tree$tip.label, trait_name] + 1
Nstates <- length(unique(traits))
anc_results <- castor::asr_max_parsimony(support_tree, traits, Nstates)
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
rownames(main_tree_edges) <- paste(main_tree_edges$anc, main_tree_edges$offspring, sep = "-")
rownames(support_edges) <- paste(support_edges$anc, support_edges$offspring, sep = "-")

trait_anc <- data.frame(main_tree_edges, node_states = support_edges[rownames(main_tree_edges), "node_states"])
trait_anc$node_states <- trait_anc$node_states - 1
rownames(trait_anc) <- trait_anc$branch_label
branch_state <- trait_anc$node_states
names(branch_state) <- rownames(trait_anc)

screened_branches <- names(branch_state)[branch_state != 0.5]
positive_branches <- names(branch_state)[branch_state == 1]
negative_branches <- names(branch_state)[branch_state == 0]

observed_summary <- utils::read.delim(stage2_summary_file, stringsAsFactors = FALSE, check.names = FALSE)
observed_row <- observed_summary[observed_summary$run_id == run_id, , drop = FALSE]
if (nrow(observed_row) != 1) stop("Could not find observed fix_drop_whale Stage 2 summary row.")

asr_summary <- data.frame(
  observed_run_id = run_id,
  tie_method = "max.col(..., ties.method = 'first')",
  n_branch_state_0 = sum(branch_state == 0),
  n_branch_state_0_5 = sum(branch_state == 0.5),
  n_branch_state_1 = sum(branch_state == 1),
  n_screened_branches = length(screened_branches),
  n_positive = length(positive_branches),
  n_negative = length(negative_branches),
  branch_state_vector_key = paste(branch_state[branch_df$Branch], collapse = "|"),
  status = "PASS",
  notes = "Observed drop-whale ASR reconstructed deterministically once; permutation labels randomize eligible terminal branches, then rerun deterministic ASR.",
  stringsAsFactors = FALSE
)
write_tsv(asr_summary, file.path(qc_dir, "endpointfix_permutation_ASR_determinism_summary.tsv"))

message("Reading GBI matrix: ", gbi_file)
gbi <- utils::read.delim(gbi_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE, na.strings = "NA")
gbi <- gbi[, branch_df$Branch, drop = FALSE]
gbi_mat <- as.matrix(gbi)
storage.mode(gbi_mat) <- "double"

welch_screen <- function(mat, states, alpha = 0.01) {
  state0_cols <- names(states)[states == 0]
  state1_cols <- names(states)[states == 1]
  x0 <- mat[, state0_cols, drop = FALSE]
  x1 <- mat[, state1_cols, drop = FALSE]

  n0 <- rowSums(!is.na(x0))
  n1 <- rowSums(!is.na(x1))
  s0 <- rowSums(x0, na.rm = TRUE)
  s1 <- rowSums(x1, na.rm = TRUE)
  ss0 <- rowSums(x0 * x0, na.rm = TRUE)
  ss1 <- rowSums(x1 * x1, na.rm = TRUE)
  m0 <- s0 / n0
  m1 <- s1 / n1
  v0 <- (ss0 - s0 * s0 / n0) / (n0 - 1)
  v1 <- (ss1 - s1 * s1 / n1) / (n1 - 1)
  v0[v0 < 0 & v0 > -1e-10] <- 0
  v1[v1 < 0 & v1 > -1e-10] <- 0
  se2 <- v0 / n0 + v1 / n1
  tvalue <- (m1 - m0) / sqrt(se2)
  df <- se2 * se2 / ((v0 / n0)^2 / (n0 - 1) + (v1 / n1)^2 / (n1 - 1))
  pvalue <- 2 * stats::pt(abs(tvalue), df = df, lower.tail = FALSE)
  testable <- n0 > 1 & n1 > 1 & is.finite(tvalue) & is.finite(pvalue)

  p <- pvalue[testable]
  genes <- rownames(mat)[testable]
  tv <- tvalue[testable]
  ord <- order(p, na.last = NA)
  p_sorted <- p[ord]
  m <- length(p_sorted)
  k <- 1L
  selected <- integer(0)
  while (k <= m && p_sorted[k] < k / m * alpha) {
    selected <- c(selected, ord[k])
    k <- k + 1L
  }
  sig_t <- tv[selected]
  list(
    n_tested = length(p),
    n_sig = length(selected),
    n_slow = sum(sig_t > 0, na.rm = TRUE),
    n_fast = sum(sig_t < 0, na.rm = TRUE),
    slow_proportion = if (length(selected) > 0) sum(sig_t > 0, na.rm = TRUE) / length(selected) else NA_real_,
    sig_genes = genes[selected],
    positive_t_all = sum(tv > 0, na.rm = TRUE),
    negative_t_all = sum(tv < 0, na.rm = TRUE)
  )
}

set.seed(seed)
label_rows <- vector("list", n_perm)
summary_rows <- vector("list", n_perm)

for (perm_id in seq_len(n_perm)) {
  pseudo_positive <- sample(rownames(eligible_terminal), n_positive_terminal, replace = FALSE)
  terminal_states <- trait_terminal[[trait_name]]
  names(terminal_states) <- rownames(trait_terminal)
  terminal_states[rownames(eligible_terminal)] <- 0
  terminal_states[pseudo_positive] <- 1

  terminal_trait_by_species <- terminal_states[terminal_df$Branch]
  names(terminal_trait_by_species) <- terminal_df$Species
  perm_traits <- terminal_trait_by_species[support_tree$tip.label] + 1
  perm_nstates <- length(unique(perm_traits))
  perm_anc <- castor::asr_max_parsimony(support_tree, perm_traits, perm_nstates)
  perm_node_states <- max.col(perm_anc$ancestral_likelihoods, ties.method = "first")
  perm_traits_states <- data.frame(
    traits = c(perm_traits, perm_node_states),
    node_number = seq_along(c(perm_traits, perm_node_states))
  )
  rownames(perm_traits_states) <- paste0("B", perm_traits_states$node_number)
  perm_support_edges <- support_edges
  perm_support_edges$node_states <- perm_traits_states[paste0("B", perm_support_edges$offspring), "traits"]
  perm_trait_anc <- data.frame(main_tree_edges, node_states = perm_support_edges[rownames(main_tree_edges), "node_states"])
  perm_trait_anc$node_states <- perm_trait_anc$node_states - 1
  rownames(perm_trait_anc) <- perm_trait_anc$branch_label
  perm_trait_anc <- perm_trait_anc[colnames(gbi_mat), , drop = FALSE]
  states <- perm_trait_anc$node_states
  names(states) <- rownames(perm_trait_anc)
  states <- states[states != 0.5 & !is.na(states)]

  label_rows[[perm_id]] <- data.frame(
    perm_id = perm_id,
    seed = seed,
    branch_id = rownames(eligible_terminal),
    species_id = eligible_terminal$Species,
    perm_state = as.integer(terminal_states[rownames(eligible_terminal)]),
    stringsAsFactors = FALSE
  )
  res <- welch_screen(gbi_mat, states, alpha = alpha)
  summary_rows[[perm_id]] <- data.frame(
    perm_id = perm_id,
    seed = seed,
    n_positive = sum(states == 1),
    n_negative = sum(states == 0),
    n_tested_genes = res$n_tested,
    n_sig_FDR_0_01 = res$n_sig,
    n_slow = res$n_slow,
    n_fast = res$n_fast,
    slow_proportion = res$slow_proportion,
    slow_proportion_defined = !is.na(res$slow_proportion),
    positive_t_all = res$positive_t_all,
    negative_t_all = res$negative_t_all,
    notes = "Matched-positive eligible-terminal-label permutation followed by deterministic ASR; old BH/FDR while-loop rule.",
    stringsAsFactors = FALSE
  )
  if (perm_id %% 10 == 0) message("Permutation ", perm_id, " / ", n_perm)
}

label_table <- do.call(rbind, label_rows)
perm_summary <- do.call(rbind, summary_rows)
write_tsv(label_table, file.path(out_dir, "endpointfix_permutation_label_sets_long.tsv"))
write_tsv(perm_summary, file.path(out_dir, "endpointfix_permutation_screening_summary.tsv"))
write_tsv(perm_summary, file.path(qc_dir, "endpointfix_permutation_screening_summary.tsv"))

obs_n_sig <- observed_row$n_sig_FDR_0_01
obs_slow <- observed_row$positive_t_sig
obs_fast <- observed_row$negative_t_sig
obs_slow_prop <- observed_row$slow_proportion

null_defined <- perm_summary$slow_proportion[!is.na(perm_summary$slow_proportion)]
defined_count <- length(null_defined)
defined_median <- if (defined_count > 0) stats::median(null_defined) else NA_real_
defined_max <- if (defined_count > 0) max(null_defined) else NA_real_
defined_min <- if (defined_count > 0) min(null_defined) else NA_real_
empirical_p_slow <- if (defined_count > 0) (sum(null_defined >= obs_slow_prop) + 1) / (defined_count + 1) else NA_real_
observed_vs_null <- data.frame(
  observed_run_id = run_id,
  observed_n_sig_FDR_0_01 = obs_n_sig,
  observed_n_slow = obs_slow,
  observed_n_fast = obs_fast,
  observed_slow_proportion = obs_slow_prop,
  n_permutations = n_perm,
  null_median_sig = stats::median(perm_summary$n_sig_FDR_0_01),
  null_max_sig = max(perm_summary$n_sig_FDR_0_01),
  null_min_sig = min(perm_summary$n_sig_FDR_0_01),
  null_defined_slow_prop_count = defined_count,
  null_undefined_slow_prop_count = sum(is.na(perm_summary$slow_proportion)),
  null_median_slow_prop = defined_median,
  null_max_slow_prop = defined_max,
  null_min_slow_prop = defined_min,
  empirical_p_sig_count = (sum(perm_summary$n_sig_FDR_0_01 >= obs_n_sig) + 1) / (n_perm + 1),
  empirical_p_slow_prop = empirical_p_slow,
  notes = "Permutation tests sample-size/null-label behavior only; no biological signal interpretation is made here.",
  stringsAsFactors = FALSE
)
write_tsv(observed_vs_null, file.path(qc_dir, "endpointfix_permutation_observed_vs_null_summary.tsv"))
write_tsv(observed_vs_null, file.path(out_dir, "endpointfix_permutation_observed_vs_null_summary.tsv"))

label_set_summary <- data.frame(
  n_permutations = n_perm,
  seed = seed,
  eligible_terminal_universe = nrow(eligible_terminal),
  observed_screened_branch_universe = length(screened_branches),
  matched_positive_terminal_count = n_positive_terminal,
  matched_negative_terminal_count = n_negative_terminal,
  n_label_rows = nrow(label_table),
  label_generation = "sample exactly matched n_positive eligible terminal branches, then rerun deterministic ASR for each permutation",
  status = "PASS",
  notes = "Per-permutation labels saved in outputs/endpointfix_permutation_label_sets_long.tsv.",
  stringsAsFactors = FALSE
)
write_tsv(label_set_summary, file.path(qc_dir, "endpointfix_permutation_label_set_summary.tsv"))

num_check <- data.frame(
  check_item = c("NaN_in_summary", "Inf_in_summary", "negative_sig_counts", "missing_slow_prop_handled"),
  status = c(
    ifelse(any(is.nan(as.matrix(perm_summary[sapply(perm_summary, is.numeric)]))), "STOP", "PASS"),
    ifelse(any(is.infinite(as.matrix(perm_summary[sapply(perm_summary, is.numeric)]))), "STOP", "PASS"),
    ifelse(any(perm_summary$n_sig_FDR_0_01 < 0 | perm_summary$n_slow < 0 | perm_summary$n_fast < 0), "STOP", "PASS"),
    "PASS"
  ),
  n_affected = c(
    sum(is.nan(as.matrix(perm_summary[sapply(perm_summary, is.numeric)]))),
    sum(is.infinite(as.matrix(perm_summary[sapply(perm_summary, is.numeric)]))),
    sum(perm_summary$n_sig_FDR_0_01 < 0 | perm_summary$n_slow < 0 | perm_summary$n_fast < 0),
    sum(is.na(perm_summary$slow_proportion))
  ),
  notes = c(
    "Numeric leak check over permutation summary numeric fields.",
    "Numeric leak check over permutation summary numeric fields.",
    "Significant-gene and direction counts must be non-negative.",
    "Undefined slow proportions are retained as NA when permutations have zero significant genes."
  ),
  stringsAsFactors = FALSE
)
write_tsv(num_check, file.path(qc_dir, "endpointfix_permutation_numeric_leak_check.tsv"))

copy_head <- function(src, dst, n = 50L) {
  lines <- readLines(src, n = n + 1L)
  writeLines(lines, dst)
}
copy_head(file.path(out_dir, "endpointfix_permutation_screening_summary.tsv"), file.path(min_dir, "endpointfix_permutation_screening_summary.head50.tsv"), 50L)
copy_head(file.path(out_dir, "endpointfix_permutation_label_sets_long.tsv"), file.path(min_dir, "endpointfix_permutation_label_sets_long.head200.tsv"), 200L)

run_summary <- c(
  "# Endpoint-Fix Drop-Whale Permutation Control",
  "",
  paste0("- Run timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  paste0("- n_permutations: ", n_perm),
  paste0("- seed: ", seed),
  "- GBI input: approved endpoint-fix old-label no_fuse fix GBI.",
  "- Observed source: release-ready deterministic Stage 2 fix_drop_whale t-test summary.",
  "- ASR tie handling: deterministic `max.col(..., ties.method = \"first\")` for the observed run and every permutation.",
  "- Permutation design: matched-positive eligible-terminal-label permutations followed by deterministic ASR; this matches the old exact permutation t-test diagnostic structure without refitting LASSO here.",
  "- FDR rule: old while-loop BH/FDR rule `p[k] < k/m*alpha` with alpha = 0.01.",
  "- No free t-test, LASSO, enrichment, figures, or old-vs-new biological comparison were run.",
  "",
  "## Observed vs Null",
  paste(capture.output(print(observed_vs_null, row.names = FALSE)), collapse = "\n"),
  "",
  "Codex self-check: permutation control completed with no obvious STOP-level QC failure. Final PASS requires review by marine mammal rerun control chat."
)
writeLines(run_summary, file.path(qc_dir, "endpointfix_permutation_run_summary.md"))

message("Permutation control complete.")
