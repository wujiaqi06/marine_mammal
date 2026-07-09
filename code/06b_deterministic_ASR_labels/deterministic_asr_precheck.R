args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[[1]] else getwd()
setwd(root_dir)

suppressPackageStartupMessages({
  library(ape)
  library(castor)
  library(digest)
})

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, quote = FALSE, sep = "\t", row.names = FALSE)
}

trait_file <- "05_trait_tables/derived/trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"
support_tree_file <- "07_ttest_screening/inputs/branch_files/mammal302.anno.BL_support.nwk"
main_tree_file <- "07_ttest_screening/inputs/branch_files/mammal302.anno.nwk"
qc_dir <- "07_ttest_screening/stage1_deterministic_asr/qc"
traits_to_check <- c("marine_binary", "aquatic_v2", "aquatic_v1", "aquatic_v3")

trait_all <- utils::read.delim(trait_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
support_tree <- ape::read.tree(support_tree_file)
main_tree <- ape::read.tree(main_tree_file)

get_asr <- function(trait_column) {
  tip_states <- trait_all[support_tree$tip.label, trait_column] + 1
  nstates <- length(unique(tip_states))
  anc_results <- castor::asr_max_parsimony(support_tree, tip_states, nstates)
  state_values <- sort(unique(trait_all[support_tree$tip.label, trait_column]))
  output_state_values <- seq_len(ncol(anc_results$ancestral_likelihoods)) - 1
  state_column_order <- data.frame(
    trait_column = trait_column,
    state_column_index = seq_len(ncol(anc_results$ancestral_likelihoods)),
    state_value = output_state_values,
    source = "Mirrors old pipeline semantics: max.col() returns a likelihood column index and the selected internal-node state is written as column_index - 1. Terminal 0.5 states are preserved separately from tip trait values.",
    status = ifelse(length(state_values) == ncol(anc_results$ancestral_likelihoods), "PASS", "STOP"),
    notes = "Required before deterministic tie handling.",
    stringsAsFactors = FALSE
  )
  list(tip_states = tip_states, anc_results = anc_results, state_values = state_values, state_column_order = state_column_order)
}

make_branch_states <- function(trait_column, asr_obj) {
  anc_lik <- asr_obj$anc_results$ancestral_likelihoods
  node_state_col <- max.col(anc_lik, ties.method = "first")
  node_states <- node_state_col - 1
  tip_states0 <- asr_obj$tip_states - 1

  traits_states <- rbind(
    data.frame(state = tip_states0, node_number = seq_along(tip_states0)),
    data.frame(state = node_states, node_number = seq(length(tip_states0) + 1, length(tip_states0) + support_tree$Nnode))
  )
  rownames(traits_states) <- paste0("B", traits_states$node_number)

  support_edges <- data.frame(support_tree$edge)
  rownames(support_edges) <- paste0("B", support_edges[, 2])
  colnames(support_edges) <- c("anc", "offspring")
  support_edges$state <- traits_states[rownames(support_edges), "state"]

  main_edges <- data.frame(main_tree$edge, seq_len(nrow(main_tree$edge)))
  colnames(main_edges) <- c("anc", "offspring", "order")
  main_edges$branch_label <- paste0("B", main_edges$order)
  rownames(support_edges) <- paste(support_edges$anc, support_edges$offspring, sep = "-")
  rownames(main_edges) <- paste(main_edges$anc, main_edges$offspring, sep = "-")

  branch_states <- support_edges[rownames(main_edges), "state"]
  names(branch_states) <- main_edges$branch_label
  branch_states
}

state_column_rows <- list()
check_rows <- list()
tied_rows <- list()

for (trait_column in traits_to_check) {
  first_asr <- get_asr(trait_column)
  state_column_rows[[trait_column]] <- first_asr$state_column_order

  tied <- apply(first_asr$anc_results$ancestral_likelihoods, 1, function(x) which(x == max(x)))
  tied_idx <- which(vapply(tied, length, integer(1)) > 1)
  if (length(tied_idx) > 0) {
    tied_rows[[trait_column]] <- do.call(rbind, lapply(tied_idx, function(i) {
      tied_cols <- tied[[i]]
      data.frame(
        trait_column = trait_column,
        node_id = paste0("N", i + length(first_asr$tip_states)),
        tied_states = paste(first_asr$state_values[tied_cols], collapse = ";"),
      selected_state = tied_cols[[1]] - 1,
      tie_method = "max.col ties.method='first'",
      state_column_order = paste(paste(seq_len(ncol(first_asr$anc_results$ancestral_likelihoods)), seq_len(ncol(first_asr$anc_results$ancestral_likelihoods)) - 1, sep = "="), collapse = ";"),
        notes = "Node id is support-tree node index. First tied likelihood column is selected deterministically.",
        stringsAsFactors = FALSE
      )
    }))
  } else {
    tied_rows[[trait_column]] <- data.frame(
      trait_column = trait_column,
      node_id = "NONE",
      tied_states = "NA",
      selected_state = "NA",
      tie_method = "max.col ties.method='first'",
      state_column_order = paste(paste(seq_len(ncol(first_asr$anc_results$ancestral_likelihoods)), seq_len(ncol(first_asr$anc_results$ancestral_likelihoods)) - 1, sep = "="), collapse = ";"),
      notes = "No tied internal nodes.",
      stringsAsFactors = FALSE
    )
  }

  for (replicate in seq_len(20)) {
    asr_obj <- get_asr(trait_column)
    branch_states <- make_branch_states(trait_column, asr_obj)
    tied_count <- sum(apply(asr_obj$anc_results$ancestral_likelihoods, 1, function(x) sum(x == max(x)) > 1))
    check_rows[[paste(trait_column, replicate, sep = "_")]] <- data.frame(
      trait_column = trait_column,
      replicate = replicate,
      tie_method = "max.col ties.method='first'",
      tied_internal_nodes = tied_count,
      n_branch_state_0 = sum(branch_states == 0),
      n_branch_state_0_5 = sum(branch_states == 0.5),
      n_branch_state_1 = sum(branch_states == 1),
      n_branch_state_other = sum(!(branch_states %in% c(0, 0.5, 1))),
      branch_state_vector_sha256 = digest::digest(paste(names(branch_states), branch_states, sep = "=", collapse = "|"), algo = "sha256", serialize = FALSE),
      status = "PENDING",
      notes = "Repeated deterministic ASR check before t-test.",
      stringsAsFactors = FALSE
    )
  }
}

asr_check <- do.call(rbind, check_rows)
asr_check$status <- ave(asr_check$branch_state_vector_sha256, asr_check$trait_column, FUN = function(x) if (length(unique(x)) == 1) "PASS" else "STOP")

summary_rows <- do.call(rbind, lapply(split(asr_check, asr_check$trait_column), function(x) {
  data.frame(
    trait_column = x$trait_column[[1]],
    n_replicates = nrow(x),
    tie_method = x$tie_method[[1]],
    tied_internal_nodes = paste(unique(x$tied_internal_nodes), collapse = ";"),
    n_unique_branch_state_sha256 = length(unique(x$branch_state_vector_sha256)),
    n_branch_state_0_range = paste(range(x$n_branch_state_0), collapse = "-"),
    n_branch_state_0_5_range = paste(range(x$n_branch_state_0_5), collapse = "-"),
    n_branch_state_1_range = paste(range(x$n_branch_state_1), collapse = "-"),
    n_branch_state_other_range = paste(range(x$n_branch_state_other), collapse = "-"),
    status = ifelse(length(unique(x$branch_state_vector_sha256)) == 1, "PASS", "STOP"),
    notes = "PASS requires identical branch state vector hashes across all 20 replicates.",
    stringsAsFactors = FALSE
  )
}))

write_tsv(do.call(rbind, state_column_rows), file.path(qc_dir, "endpointfix_stage1_asr_state_column_order.tsv"))
write_tsv(asr_check, file.path(qc_dir, "endpointfix_stage1_deterministic_asr_check.tsv"))
write_tsv(summary_rows, file.path(qc_dir, "endpointfix_stage1_deterministic_asr_summary.tsv"))
write_tsv(do.call(rbind, tied_rows), file.path(qc_dir, "endpointfix_stage1_tied_internal_nodes.tsv"))

if (any(summary_rows$status != "PASS") || any(do.call(rbind, state_column_rows)$status != "PASS")) {
  stop("Deterministic ASR precheck failed.")
}
