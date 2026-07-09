args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript run_ttest_screening.R <manifest.tsv> <run_id> <out_dir>")
}

manifest_file <- args[1]
run_id <- args[2]
out_dir <- args[3]

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

script_path_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- sub("^--file=", "", script_path_arg[1] %||% "run_ttest_screening.R")
script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
source(file.path(script_dir, "marine_functions.R"))

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

manifest <- utils::read.delim(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!"run_id" %in% colnames(manifest)) {
  stop("Manifest must contain run_id.")
}
run_row <- manifest[manifest$run_id == run_id, , drop = FALSE]
if (nrow(run_row) != 1) {
  stop(sprintf("Expected exactly one manifest row for run_id '%s'.", run_id))
}

required_fields <- c(
  "trait_name",
  "trait_file",
  "trait_column",
  "gbi_matrix",
  "branch_file",
  "tree_file",
  "support_tree_file",
  "alpha",
  "out_prefix"
)
missing_fields <- setdiff(required_fields, colnames(run_row))
if (length(missing_fields) > 0) {
  stop(sprintf("Manifest is missing required fields: %s", paste(missing_fields, collapse = ", ")))
}

trait_name <- run_row$trait_name[[1]]
trait_file <- run_row$trait_file[[1]]
trait_column <- run_row$trait_column[[1]]
gbi_matrix_file <- run_row$gbi_matrix[[1]]
branch_file <- run_row$branch_file[[1]]
tree_file <- run_row$tree_file[[1]]
support_tree_file <- run_row$support_tree_file[[1]]
out_prefix <- run_row$out_prefix[[1]]
alpha <- as.numeric(run_row$alpha[[1]])
clade_name <- "mammal"

trait_all <- utils::read.delim(trait_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
if (!trait_column %in% colnames(trait_all)) {
  stop(sprintf("Trait column '%s' not found in trait file.", trait_column))
}
if (trait_name != trait_column) {
  stop(sprintf("Manifest mismatch: trait_name '%s' does not match trait_column '%s'.", trait_name, trait_column))
}

branch_df <- utils::read.delim(branch_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Species", "Branch") %in% colnames(branch_df))) {
  stop("Branch file must contain Species and Branch columns.")
}
assert_unique_ids(branch_df$Branch, "branch IDs in branch file")
branch_order <- branch_df$Branch
terminal_df <- branch_df[branch_df$Species != "internal", , drop = FALSE]

support_tree <- ape::read.tree(support_tree_file)
main_tree <- ape::read.tree(tree_file)
missing_tips <- setdiff(support_tree$tip.label, rownames(trait_all))
if (length(missing_tips) > 0) {
  stop(sprintf("Trait file is missing tip labels from support tree: %s", paste(head(missing_tips, 10), collapse = ", ")))
}

if (!file.exists(gbi_matrix_file)) {
  stop(sprintf("GBI matrix file not found: %s", gbi_matrix_file))
}
gene_branch_interaction <- utils::read.table(gbi_matrix_file, row.names = 1, check.names = FALSE, header = TRUE)
if (is.null(rownames(gene_branch_interaction)) || any(rownames(gene_branch_interaction) == "")) {
  stop("GBI matrix must have non-empty gene rownames.")
}
assert_unique_ids(rownames(gene_branch_interaction), "gene IDs in GBI matrix")
assert_unique_ids(colnames(gene_branch_interaction), "branch IDs in GBI matrix")

missing_branches <- setdiff(branch_order, colnames(gene_branch_interaction))
if (length(missing_branches) > 0) {
  stop(sprintf("Branch file contains branches not present in GBI matrix: %s", paste(head(missing_branches, 10), collapse = ", ")))
}
gene_branch_interaction <- gene_branch_interaction[, branch_order, drop = FALSE]

traits <- trait_all[support_tree$tip.label, trait_column] + 1
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
ape_bl_match <- main_tree_edges
rownames(ape_bl_match) <- paste(ape_bl_match$anc, ape_bl_match$offspring, sep = "-")
utils::write.table(
  ape_bl_match,
  file = file.path(out_dir, paste0(out_prefix, ".", clade_name, ".Ape_Bl_march.txt")),
  quote = FALSE,
  row.names = FALSE,
  sep = "\t"
)

rownames(support_edges) <- paste(support_edges$anc, support_edges$offspring, sep = "-")
trait_anc <- data.frame(ape_bl_match, node_states = support_edges[rownames(ape_bl_match), 3])
trait_anc$node_states <- trait_anc$node_states - 1
rownames(trait_anc) <- trait_anc$branch_label
trait_anc_screen <- trait_anc[trait_anc$node_states != 0.5, , drop = FALSE]
trait_anc_screen <- trait_anc_screen[rownames(trait_anc_screen) %in% colnames(gene_branch_interaction), , drop = FALSE]

branch_trait <- trait_anc[colnames(gene_branch_interaction), "node_states"]
gene_branch_interaction_trait <- rbind(trait = branch_trait, gene_branch_interaction)
gene_branch_interaction_trait <- gene_branch_interaction_trait[, rownames(trait_anc_screen), drop = FALSE]

trait_t <- list()
out_idx <- 0L
for (i in 2:nrow(gene_branch_interaction_trait)) {
  data1 <- data.frame(gene_branch_interaction_trait[c(1, i), ], check.names = FALSE)
  data11 <- t(stats::na.omit(t(data1)))
  trait_count <- table(data11[1, ])
  if ((length(trait_count) > 1) && (min(trait_count) > 1)) {
    trait_t0 <- stats::t.test(data11[2, ] ~ data11[1, ])
    out_idx <- out_idx + 1L
    trait_t[[out_idx]] <- data.frame(
      gene = rownames(gene_branch_interaction_trait)[i],
      tvalue = unname(trait_t0$statistic),
      pvalue = trait_t0$p.value,
      mean1 = unname(trait_t0$estimate[1]),
      mean2 = unname(trait_t0$estimate[2]),
      marine = unname(trait_count[2]),
      non_marine = unname(trait_count[1]),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }
}
trait_t <- do.call(rbind, trait_t)
rownames(trait_t) <- trait_t$gene
trait_t <- trait_t[order(trait_t$pvalue, decreasing = FALSE), , drop = FALSE]

ttest_file <- file.path(out_dir, paste0(out_prefix, ".t_test.n", nrow(trait_t), ".txt"))
utils::write.table(trait_t, file = ttest_file, quote = FALSE, row.names = FALSE, sep = "\t")

trait_pvalue_sort <- sort(trait_t$pvalue, na.last = NA)
trait_pvalue_order <- order(trait_t$pvalue, na.last = NA)
gene_list <- integer(0)
m <- length(trait_pvalue_sort)
if (m > 0) {
  k <- 1L
  while (k <= m && trait_pvalue_sort[k] < k / m * alpha) {
    gene_list <- c(gene_list, trait_pvalue_order[k])
    k <- k + 1L
  }
}

t_gene <- trait_t[gene_list, , drop = FALSE]
fdr_file <- file.path(out_dir, paste0(out_prefix, ".", clade_name, ".FDR", alpha, ".n", length(gene_list), ".t_test.txt"))
utils::write.table(t_gene, file = fdr_file, quote = FALSE, sep = "\t", row.names = FALSE)

t_gene_sig <- trait_t[t_gene$gene, , drop = FALSE]
t_gene_sig <- t_gene_sig[order(t_gene_sig$tvalue), , drop = FALSE]

gene_branch_interaction_sig <- t(gene_branch_interaction[t_gene_sig$gene, , drop = FALSE])
gene_branch_interaction_sig_terminal <- as.matrix(gene_branch_interaction_sig[terminal_df$Branch, , drop = FALSE])
terminal_mean <- apply(gene_branch_interaction_sig_terminal, 2, mean_no_na)
rate_all_imputed <- impute_by_column_terminal_mean(gene_branch_interaction_sig, terminal_mean)
imputed_file <- file.path(out_dir, paste0(out_prefix, ".FDR", alpha, ".GBI.imputated.n", length(gene_list), ".txt"))
utils::write.table(rate_all_imputed, file = imputed_file, quote = FALSE, sep = "\t")

run_metadata <- data.frame(
  field = c(
    "run_id",
    "trait_name",
    "trait_file",
    "trait_column",
    "gbi_matrix",
    "branch_file",
    "tree_file",
    "support_tree_file",
    "alpha",
    "out_dir",
    "timestamp"
  ),
  value = c(
    run_id,
    trait_name,
    normalizePath(trait_file, mustWork = FALSE),
    trait_column,
    normalizePath(gbi_matrix_file, mustWork = FALSE),
    normalizePath(branch_file, mustWork = FALSE),
    normalizePath(tree_file, mustWork = FALSE),
    normalizePath(support_tree_file, mustWork = FALSE),
    alpha,
    normalizePath(out_dir, mustWork = FALSE),
    format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
  ),
  stringsAsFactors = FALSE
)
utils::write.table(run_metadata, file = file.path(out_dir, "run_metadata.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

screening_summary <- data.frame(
  run_id = run_id,
  trait_name = trait_name,
  n_branches_input = ncol(gene_branch_interaction),
  n_branches_screened = ncol(gene_branch_interaction_trait),
  n_tested_genes = nrow(trait_t),
  n_sig_fdr01 = length(gene_list),
  positive_t_all = sum(trait_t$tvalue > 0, na.rm = TRUE),
  negative_t_all = sum(trait_t$tvalue < 0, na.rm = TRUE),
  positive_t_sig = sum(t_gene_sig$tvalue > 0, na.rm = TRUE),
  negative_t_sig = sum(t_gene_sig$tvalue < 0, na.rm = TRUE),
  stringsAsFactors = FALSE
)
utils::write.table(screening_summary, file = file.path(out_dir, "screening_summary.tsv"), quote = FALSE, sep = "\t", row.names = FALSE)

message("Screening complete: ", out_dir)
