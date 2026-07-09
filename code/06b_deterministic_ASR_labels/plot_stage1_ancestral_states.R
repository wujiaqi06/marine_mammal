args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) args[[1]] else getwd()
setwd(root_dir)

suppressPackageStartupMessages({
  library(ape)
  library(castor)
})

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, quote = FALSE, sep = "\t", row.names = FALSE)
}

trait_file <- "05_trait_tables/derived/trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"
support_tree_file <- "07_ttest_screening/inputs/branch_files/mammal302.anno.BL_support.nwk"
main_tree_file <- "07_ttest_screening/inputs/branch_files/mammal302.anno.nwk"
branch_file <- "07_ttest_screening/inputs/branch_files/mammal.branch.txt"
out_dir <- "07_ttest_screening/stage1_deterministic_asr/ancestral_state_plots"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

runs <- data.frame(
  run_id = c("fix_marine_binary", "fix_aquatic_v2", "fix_aquatic_v1", "fix_aquatic_v3"),
  trait_column = c("marine_binary", "aquatic_v2", "aquatic_v1", "aquatic_v3"),
  plot_title = c(
    "fix_marine_binary deterministic ancestral states",
    "fix_aquatic_v2 deterministic ancestral states",
    "fix_aquatic_v1 deterministic ancestral states",
    "fix_aquatic_v3 deterministic ancestral states"
  ),
  stringsAsFactors = FALSE
)

trait_all <- utils::read.delim(trait_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
support_tree <- ape::read.tree(support_tree_file)
main_tree <- ape::read.tree(main_tree_file)
branch_df <- utils::read.delim(branch_file, stringsAsFactors = FALSE, check.names = FALSE)

if (!all(support_tree$tip.label %in% rownames(trait_all))) {
  stop("Trait table is missing support-tree tips.")
}

make_main_edge_map <- function() {
  main_edges <- data.frame(
    anc = main_tree$edge[, 1],
    offspring = main_tree$edge[, 2],
    branch_label_raw = main_tree$edge.length,
    order = seq_len(nrow(main_tree$edge)),
    stringsAsFactors = FALSE
  )
  main_edges$branch_label <- paste0("B", as.integer(round(main_edges$branch_label_raw)))
  main_edges$row_key <- paste(main_edges$anc, main_edges$offspring, sep = "-")
  main_edges
}

make_support_edge_map <- function() {
  support_edges <- data.frame(
    anc = support_tree$edge[, 1],
    offspring = support_tree$edge[, 2],
    support_edge_order = seq_len(nrow(support_tree$edge)),
    stringsAsFactors = FALSE
  )
  support_edges$row_key <- paste(support_edges$anc, support_edges$offspring, sep = "-")
  support_edges
}

main_edges <- make_main_edge_map()
support_edges_template <- make_support_edge_map()
if (!all(support_edges_template$row_key %in% main_edges$row_key)) {
  stop("Support-tree edges cannot be matched to old-label main-tree edges by node pair.")
}
branch_label_by_row_key <- setNames(main_edges$branch_label, main_edges$row_key)

state_color <- function(x) {
  ifelse(
    is.na(x), "#BDBDBD",
    ifelse(
      x == 0, "#D9D9D9",
      ifelse(x == 0.5, "#F1A340",
        ifelse(x == 1, "#2C7FB8", "#7B3294")
      )
    )
  )
}

state_label <- function(x) {
  ifelse(
    is.na(x), "NA",
    ifelse(
      x == 0, "0",
      ifelse(x == 0.5, "0.5", ifelse(x == 1, "1", paste0("other:", x)))
    )
  )
}

get_asr <- function(trait_column) {
  # This mirrors the deterministic Stage 1 t-test ASR layer.
  tip_states <- trait_all[support_tree$tip.label, trait_column] + 1
  nstates <- length(unique(tip_states))
  anc_results <- castor::asr_max_parsimony(support_tree, tip_states, nstates)
  node_state_col <- max.col(anc_results$ancestral_likelihoods, ties.method = "first")
  node_states <- node_state_col - 1
  tip_states0 <- tip_states - 1

  all_states <- rbind(
    data.frame(state = tip_states0, node_number = seq_along(tip_states0), node_type = "terminal"),
    data.frame(
      state = node_states,
      node_number = seq(length(tip_states0) + 1, length(tip_states0) + support_tree$Nnode),
      node_type = "internal"
    )
  )
  rownames(all_states) <- paste0("B", all_states$node_number)

  support_edges <- support_edges_template
  support_edges$support_child_node_id <- paste0("B", support_edges$offspring)
  support_edges$old_branch_label <- branch_label_by_row_key[support_edges$row_key]
  support_edges$state <- all_states[support_edges$support_child_node_id, "state"]
  support_edges$node_type <- all_states[support_edges$support_child_node_id, "node_type"]
  support_edges
}

plot_one_trait <- function(run_id, trait_column, plot_title, pdf_path) {
  asr_edges <- get_asr(trait_column)
  edge_colors <- state_color(asr_edges$state)

  tip_trait <- trait_all[support_tree$tip.label, trait_column]
  tip_colors <- state_color(tip_trait)

  grDevices::pdf(pdf_path, width = 8.5, height = 11, useDingbats = FALSE)
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit({
    graphics::par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  graphics::par(mar = c(0.5, 0.5, 3.2, 0.5), xpd = NA)
  ape::plot.phylo(
    support_tree,
    no.margin = TRUE,
    edge.color = edge_colors,
    direction = "r",
    align.tip.label = FALSE,
    cex = 0.18,
    tip.color = tip_colors
  )
  graphics::title(main = plot_title, cex.main = 1.1, line = 1.4)
  graphics::legend(
    "topright",
    legend = c("state 0", "state 0.5", "state 1", "other/unexpected"),
    col = c("#D9D9D9", "#F1A340", "#2C7FB8", "#7B3294"),
    lwd = 5,
    bty = "n",
    cex = 0.75,
    title = "Branch / tip state"
  )
  graphics::mtext(
    "Deterministic ASR: castor::asr_max_parsimony + max.col(..., ties.method = 'first')",
    side = 1,
    adj = 0,
    line = -1,
    cex = 0.55
  )
  invisible(asr_edges)
}

all_branch_tables <- list()
summary_rows <- list()
for (i in seq_len(nrow(runs))) {
  run_id <- runs$run_id[[i]]
  trait_column <- runs$trait_column[[i]]
  pdf_path <- file.path(out_dir, paste0(run_id, ".ancestral_state.pdf"))
  asr_edges <- plot_one_trait(run_id, trait_column, runs$plot_title[[i]], pdf_path)

  branch_table <- merge(
    asr_edges[, c("support_edge_order", "old_branch_label", "support_child_node_id", "state", "node_type")],
    branch_df,
    by.x = "old_branch_label",
    by.y = "Branch",
    all.x = TRUE,
    sort = FALSE
  )
  branch_table$run_id <- run_id
  branch_table$trait_column <- trait_column
  branch_table$state_label <- state_label(branch_table$state)
  branch_table$plot_pdf <- normalizePath(pdf_path, mustWork = FALSE)
  branch_table <- branch_table[, c(
    "run_id", "trait_column", "old_branch_label", "Species", "node_type",
    "support_child_node_id", "support_edge_order", "state", "state_label", "plot_pdf"
  )]
  all_branch_tables[[run_id]] <- branch_table

  summary_rows[[run_id]] <- data.frame(
    run_id = run_id,
    trait_column = trait_column,
    n_branches = nrow(branch_table),
    n_state_0 = sum(branch_table$state == 0, na.rm = TRUE),
    n_state_0_5 = sum(branch_table$state == 0.5, na.rm = TRUE),
    n_state_1 = sum(branch_table$state == 1, na.rm = TRUE),
    n_state_other = sum(!(branch_table$state %in% c(0, 0.5, 1)), na.rm = TRUE),
    pdf_file = normalizePath(pdf_path, mustWork = FALSE),
    stringsAsFactors = FALSE
  )
}

branch_states_all <- do.call(rbind, all_branch_tables)
summary_all <- do.call(rbind, summary_rows)

write_tsv(branch_states_all, file.path(out_dir, "endpointfix_stage1_ancestral_state_branch_table.tsv"))
write_tsv(summary_all, file.path(out_dir, "endpointfix_stage1_ancestral_state_plot_summary.tsv"))

combined_pdf <- file.path(out_dir, "endpointfix_stage1_ancestral_state_four_traits.pdf")
grDevices::pdf(combined_pdf, width = 8.5, height = 11, useDingbats = FALSE)
for (i in seq_len(nrow(runs))) {
  asr_edges <- get_asr(runs$trait_column[[i]])
  tip_trait <- trait_all[support_tree$tip.label, runs$trait_column[[i]]]
  graphics::par(mar = c(0.5, 0.5, 3.2, 0.5), xpd = NA)
  ape::plot.phylo(
    support_tree,
    no.margin = TRUE,
    edge.color = state_color(asr_edges$state),
    direction = "r",
    align.tip.label = FALSE,
    cex = 0.18,
    tip.color = state_color(tip_trait)
  )
  graphics::title(main = runs$plot_title[[i]], cex.main = 1.1, line = 1.4)
  graphics::legend(
    "topright",
    legend = c("state 0", "state 0.5", "state 1", "other/unexpected"),
    col = c("#D9D9D9", "#F1A340", "#2C7FB8", "#7B3294"),
    lwd = 5,
    bty = "n",
    cex = 0.75,
    title = "Branch / tip state"
  )
}
grDevices::dev.off()

run_summary <- file.path(out_dir, "endpointfix_stage1_ancestral_state_plot_run_summary.md")
writeLines(c(
  "# Stage 1 Deterministic Ancestral-State Plot Summary",
  "",
  "Scope: diagnostic ancestral-state plots for the four deterministic Stage 1 baseline/comparison traits only.",
  "",
  "No t-test, LASSO, enrichment, figure assembly, or old-vs-new biological comparison was run.",
  "",
  paste0("Trait table: `", normalizePath(trait_file, mustWork = FALSE), "`"),
  paste0("Support tree: `", normalizePath(support_tree_file, mustWork = FALSE), "`"),
  paste0("Old-label tree: `", normalizePath(main_tree_file, mustWork = FALSE), "`"),
  paste0("Branch table: `", normalizePath(branch_file, mustWork = FALSE), "`"),
  "",
  "ASR logic: `castor::asr_max_parsimony()` followed by `max.col(..., ties.method = \"first\")`, matching deterministic Stage 1.",
  "",
  "Generated PDFs:",
  paste0("- `", summary_all$pdf_file, "`"),
  paste0("- `", normalizePath(combined_pdf, mustWork = FALSE), "`"),
  "",
  "Branch-state summary:",
  paste(apply(summary_all[, c("run_id", "n_state_0", "n_state_0_5", "n_state_1", "n_state_other")], 1, paste, collapse = "\t"), collapse = "\n"),
  "",
  "Codex self-check: ancestral-state diagnostic PDFs generated from deterministic Stage 1 ASR logic."
), con = run_summary)

message("Ancestral-state plots written to: ", normalizePath(out_dir, mustWork = FALSE))
