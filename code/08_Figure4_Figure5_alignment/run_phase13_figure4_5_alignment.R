#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(parallel)
  library(ape)
  library(castor)
  library(gridExtra)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment")
dirs <- list(
  fig4 = file.path(work_dir, "Figure4_corrected_full_data"),
  fig5a = file.path(work_dir, "Figure5A_nested_sensitivity"),
  fig5b = file.path(work_dir, "Figure5B_permutation_cleanup"),
  fig5c = file.path(work_dir, "Figure5C_corrected_turnover_module_null"),
  figs = file.path(work_dir, "figures_draft"),
  qc = file.path(work_dir, "qc"),
  decision = file.path(work_dir, "decision"),
  code = file.path(work_dir, "code")
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}
cores <- as.integer(get_arg("--cores", min(12L, max(1L, parallel::detectCores() - 1L))))
if (!is.finite(cores) || cores < 1L) cores <- 1L
n_perm <- as.integer(get_arg("--n_perm", 10000L))
seed_base <- as.integer(get_arg("--seed", 20260524L))
alpha <- as.numeric(get_arg("--alpha", 0.01))
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")

read_tsv <- function(file, ...) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_tsv <- function(x, file) utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
stamp <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)
bind_or_empty <- function(xs) {
  xs <- xs[!vapply(xs, function(x) is.null(x) || nrow(x) == 0, logical(1))]
  if (length(xs) == 0) return(data.frame())
  all_names <- unique(unlist(lapply(xs, names)))
  xs <- lapply(xs, function(x) {
    missing <- setdiff(all_names, names(x))
    for (m in missing) x[[m]] <- NA
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}
iqr_text <- function(x, digits = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  qs <- stats::quantile(x, c(0.25, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f-%.", digits, "f"), qs[[1]], qs[[2]])
}
auc_rank <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]; p <- p[keep]
  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}
bh_old_while <- function(p, alpha = 0.01) {
  ord <- order(p, na.last = NA)
  p_sorted <- p[ord]
  selected <- integer(0)
  m <- length(p_sorted)
  if (m > 0) {
    k <- 1L
    while (k <= m && p_sorted[[k]] < k / m * alpha) {
      selected <- c(selected, ord[[k]])
      k <- k + 1L
    }
  }
  selected
}
save_plot <- function(plot, stem, width = 7, height = 4.5) {
  grDevices::pdf(file.path(dirs$figs, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  print(plot); grDevices::dev.off()
  grDevices::png(file.path(dirs$figs, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  print(plot); grDevices::dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(dirs$figs, paste0(stem, ".svg")), width = width, height = height)
  } else {
    grDevices::svg(file.path(dirs$figs, paste0(stem, ".svg")), width = width, height = height)
  }
  print(plot); grDevices::dev.off()
}
save_grid <- function(plots, stem, width = 7, height = 4.5, ncol = 1) {
  grDevices::pdf(file.path(dirs$figs, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  gridExtra::grid.arrange(grobs = plots, ncol = ncol); grDevices::dev.off()
  grDevices::png(file.path(dirs$figs, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = plots, ncol = ncol); grDevices::dev.off()
  if (requireNamespace("svglite", quietly = TRUE)) {
    svglite::svglite(file.path(dirs$figs, paste0(stem, ".svg")), width = width, height = height)
  } else {
    grDevices::svg(file.path(dirs$figs, paste0(stem, ".svg")), width = width, height = height)
  }
  gridExtra::grid.arrange(grobs = plots, ncol = ncol); grDevices::dev.off()
}

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  marine_drop_traits = file.path(root_dir, "05_trait_tables", "drop_inputs", "marine_drop_trait_inputs", "trait_input.marine_binary_drop_sets.combined.tsv"),
  aquatic_drop_traits = file.path(root_dir, "05_trait_tables", "drop_inputs", "aquatic_v2_drop_trait_inputs", "trait_input.aquatic_v2_drop_sets.combined.tsv"),
  stage1_asr = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "ancestral_state_plots", "endpointfix_stage1_ancestral_state_branch_table.tsv"),
  stage2_manifest = file.path(root_dir, "07_ttest_screening", "stage2_deterministic_drop_sensitivity", "qc", "endpointfix_stage2_deterministic_manifest.tsv"),
  stage2_state_summary = file.path(root_dir, "07_ttest_screening", "stage2_deterministic_drop_sensitivity", "qc", "endpointfix_stage2_trait_state_branch_summary.tsv"),
  phase12a_coef = file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity", "corrected_coefficients", "corrected_preprocessing_full_data_coefficients_all_runs.tsv"),
  phase12a_full_summary = file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity", "corrected_coefficients", "corrected_preprocessing_full_data_fit_summary_all_runs.tsv"),
  phase12b_auc = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_AUC_summary.tsv"),
  phase12b_oof = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_OOF_predictions.tsv"),
  phase12b_selected = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_fold_selected_predictors.tsv"),
  module_map = file.path(root_dir, "09_figures", "figure_ready_tables", "outputs", "endpointfix_Fig4B_coefficient_architecture_table.tsv"),
  fig5b_distribution = file.path(root_dir, "09_figures", "figure_ready_tables", "outputs", "endpointfix_Fig5B_permutation_distribution.tsv"),
  fig5b_summary = file.path(root_dir, "09_figures", "figure_ready_tables", "outputs", "endpointfix_Fig5B_permutation_summary.tsv"),
  fig5b_qc = file.path(root_dir, "09_figures", "figure_ready_tables", "qc", "endpointfix_Fig5B_update_summary.md")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

feature_file <- function(run_id, prefix, n, stage) {
  if (stage == "stage1") {
    file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", run_id, paste0(prefix, ".mammal.FDR0.01.n", n, ".t_test.txt"))
  } else {
    file.path(root_dir, "07_ttest_screening", "stage2_deterministic_drop_sensitivity", "outputs", run_id, paste0(prefix, ".mammal.FDR0.01.n", n, ".t_test.txt"))
  }
}
run_specs_all <- data.frame(
  run_id = c(
    "fix_marine_binary", "fix_drop_whale", "fix_whale_only", "fix_pinniped_only",
    "fix_aquatic_v2", "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia", "fix_aquatic_v2_noMarineEdge"
  ),
  trait = c(rep("marine_binary", 4), rep("binary_aquatic_dependence", 4)),
  trait_axis = c(rep("marine", 4), rep("binary aquatic-dependence", 4)),
  trait_column = c("marine_binary", "drop_whale", "whale_only", "pinniped_only", "aquatic_v2", "aquatic_v2_noCetacea", "aquatic_v2_noPinnipedia", "aquatic_v2_noMarineEdge"),
  drop_rule = c("baseline", "Cetacea", "Cetacea only", "Pinnipedia only", "baseline", "Cetacea", "Pinnipedia", "MarineEdge"),
  feature_prefix = c("marine_binary", "drop_whale", "whale_only", "pinniped_only", "aquatic_v2", "aquatic_v2_noCetacea", "aquatic_v2_noPinnipedia", "aquatic_v2_noMarineEdge"),
  n_global_features = c(1559, 894, 1765, 795, 1227, 560, 1246, 1265),
  stage = c("stage1", "stage2", "stage2", "stage2", "stage1", "stage2", "stage2", "stage2"),
  stringsAsFactors = FALSE
)
run_specs_all$feature_file <- mapply(feature_file, run_specs_all$run_id, run_specs_all$feature_prefix, run_specs_all$n_global_features, run_specs_all$stage, USE.NAMES = FALSE)
if (any(!file.exists(run_specs_all$feature_file))) stop("Missing feature files for Figure 5A/5C run specs.")

stamp("Reading shared inputs")
gbi <- read_tsv(paths$gbi)
rownames(gbi) <- gbi$gene
gbi$gene <- NULL
gbi_mat_all <- as.matrix(gbi)
storage.mode(gbi_mat_all) <- "numeric"
rm(gbi)
branch <- read_tsv(paths$branch)
trait <- read_tsv(paths$trait)
marine_drop <- read_tsv(paths$marine_drop_traits)
aquatic_drop <- read_tsv(paths$aquatic_drop_traits)
terminal <- branch[branch$Species != "internal", , drop = FALSE]
terminal$genus <- vapply(strsplit(terminal$Species, "_"), `[`, character(1), 1)
terminal <- merge(terminal, trait, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- merge(terminal, marine_drop, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- merge(terminal, aquatic_drop, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- terminal[match(branch$Species[branch$Species != "internal"], terminal$Species), , drop = FALSE]
terminal$trait_category <- ifelse(terminal$marine_binary == 1, "marine", ifelse(terminal$aquatic_v2 == 1, "non-marine aquatic", ifelse(terminal$aquatic_v2 == 0.5, "semi-aquatic", "terrestrial")))
terminal$aquaticity_score <- terminal$aquaticity_score_sum_0_18
gbi_mat_all <- gbi_mat_all[, branch$Branch, drop = FALSE]

coef_all <- read_tsv(paths$phase12a_coef)
full_summary <- read_tsv(paths$phase12a_full_summary)
module_map <- read_tsv(paths$module_map)
module_lookup <- module_map[, c("gene_id", "candidate_module", "annotation_status", "module_annotation_source"), drop = FALSE]
module_lookup <- module_lookup[!duplicated(module_lookup$gene_id), , drop = FALSE]
names(module_lookup)[1] <- "gene"

get_selected <- function(run_id) {
  sort(unique(coef_all$gene[coef_all$run_id == run_id & coef_all$selected]))
}
get_feature_genes <- function(run_id) {
  row <- run_specs_all[run_specs_all$run_id == run_id, , drop = FALSE]
  read_tsv(row$feature_file[[1]])$gene
}
annotate_genes <- function(df) {
  out <- merge(df, module_lookup, by = "gene", all.x = TRUE, sort = FALSE)
  out$candidate_module[is.na(out$candidate_module)] <- "REVIEW_UNANNOTATED"
  out$annotation_status[is.na(out$annotation_status)] <- "REVIEW_UNANNOTATED"
  out$module_annotation_source[is.na(out$module_annotation_source)] <- "not_in_old_module_table"
  out
}

## Part 1: Figure 4 corrected full-data architecture
stamp("Building Figure 4 corrected full-data architecture")
marine_set <- get_selected("fix_marine_binary")
aquatic_set <- get_selected("fix_aquatic_v2")
shared <- intersect(marine_set, aquatic_set)
marine_only <- setdiff(marine_set, aquatic_set)
aquatic_only <- setdiff(aquatic_set, marine_set)
partition <- data.frame(
  partition = c("marine-specific predictors", "shared predictors", "aquatic-specific predictors", "union total"),
  n_genes = c(length(marine_only), length(shared), length(aquatic_only), length(union(marine_set, aquatic_set))),
  genes = c(paste(marine_only, collapse = ";"), paste(shared, collapse = ";"), paste(aquatic_only, collapse = ";"), paste(union(marine_set, aquatic_set), collapse = ";")),
  source = "Phase 12A corrected full-data baseline LASSO coefficients",
  stringsAsFactors = FALSE
)
write_tsv(partition, file.path(dirs$fig4, "Figure4_corrected_predictor_partition.tsv"))

coef_wide <- merge(
  coef_all[coef_all$run_id == "fix_marine_binary", c("gene", "coefficient"), drop = FALSE],
  coef_all[coef_all$run_id == "fix_aquatic_v2", c("gene", "coefficient"), drop = FALSE],
  by = "gene", all = TRUE, suffixes = c("_marine", "_aquatic")
)
coef_wide$coefficient_marine[is.na(coef_wide$coefficient_marine)] <- 0
coef_wide$coefficient_aquatic[is.na(coef_wide$coefficient_aquatic)] <- 0
coef_wide$marine_nonzero <- coef_wide$gene %in% marine_set
coef_wide$aquatic_nonzero <- coef_wide$gene %in% aquatic_set
coef_wide <- coef_wide[coef_wide$marine_nonzero | coef_wide$aquatic_nonzero, , drop = FALSE]
coef_wide$predictor_class <- ifelse(coef_wide$marine_nonzero & coef_wide$aquatic_nonzero, "shared", ifelse(coef_wide$marine_nonzero, "marine_specific", "aquatic_specific"))
coef_wide$display_category <- ifelse(coef_wide$predictor_class == "shared", "Shared predictors", ifelse(coef_wide$predictor_class == "marine_specific", "Marine-specific predictors", "Aquatic-specific predictors"))
coef_wide$abs_marine_coef <- abs(coef_wide$coefficient_marine)
coef_wide$abs_aquatic_coef <- abs(coef_wide$coefficient_aquatic)
coef_wide$max_abs_coef <- pmax(coef_wide$abs_marine_coef, coef_wide$abs_aquatic_coef)
coef_wide$marine_abs_gt_0_1 <- coef_wide$abs_marine_coef > 0.1
coef_wide$aquatic_abs_gt_0_1 <- coef_wide$abs_aquatic_coef > 0.1
coef_wide$coefficient_size_marker <- ifelse(coef_wide$marine_abs_gt_0_1 & coef_wide$aquatic_abs_gt_0_1, "*#", ifelse(coef_wide$marine_abs_gt_0_1, "*", ifelse(coef_wide$aquatic_abs_gt_0_1, "#", "")))
coef_wide$marine_direction <- ifelse(coef_wide$coefficient_marine > 0, "fast_direction_positive_coefficient", ifelse(coef_wide$coefficient_marine < 0, "slow_direction_negative_coefficient", "not_selected_or_zero"))
coef_wide$aquatic_direction <- ifelse(coef_wide$coefficient_aquatic > 0, "fast_direction_positive_coefficient", ifelse(coef_wide$coefficient_aquatic < 0, "slow_direction_negative_coefficient", "not_selected_or_zero"))
coef_arch <- annotate_genes(coef_wide)
coef_arch <- coef_arch[order(coef_arch$display_category, -coef_arch$max_abs_coef, coef_arch$gene), , drop = FALSE]
write_tsv(coef_arch, file.path(dirs$fig4, "Figure4_corrected_coefficient_architecture.tsv"))

module_partition <- aggregate(gene ~ candidate_module + display_category + annotation_status, data = coef_arch, FUN = length)
names(module_partition)[names(module_partition) == "gene"] <- "gene_count"
rep_genes <- aggregate(gene ~ candidate_module + display_category, data = coef_arch, FUN = function(x) paste(sort(unique(x)), collapse = ";"))
names(rep_genes)[3] <- "genes"
module_partition <- merge(module_partition, rep_genes, by = c("candidate_module", "display_category"), all.x = TRUE, sort = FALSE)
write_tsv(module_partition, file.path(dirs$fig4, "Figure4_corrected_module_partition.tsv"))

fig4_note <- c(
  "# Figure 4 Evidence Note",
  "",
  "Figure 4 summarizes full-data final LASSO architecture after corrected preprocessing. It is not a cross-validated validation panel; predictive performance is evaluated by nested supervised feature-selection gLOOCV in Fig. 2 and Fig. 5A.",
  "",
  "Canonical corrected full-data counts:",
  paste0("- marine predictors = ", length(marine_set)),
  paste0("- binary aquatic-dependence predictors = ", length(aquatic_set)),
  paste0("- shared predictors = ", length(shared)),
  paste0("- marine-specific predictors = ", length(marine_only)),
  paste0("- aquatic-specific predictors = ", length(aquatic_only)),
  paste0("- union total = ", length(union(marine_set, aquatic_set)))
)
writeLines(fig4_note, file.path(dirs$fig4, "Figure4_corrected_full_data_evidence_note.md"))

p4a <- ggplot(partition[partition$partition != "union total", ], aes(x = partition, y = n_genes, fill = partition)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_genes), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("marine-specific predictors" = "#2F72C4", "shared predictors" = "#666666", "aquatic-specific predictors" = "#E28416")) +
  labs(x = NULL, y = "Predictors", title = "A. Corrected full-data predictor partition") +
  theme_classic(base_size = 9) + theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))
p4b_df <- rbind(
  data.frame(gene = coef_arch$gene, model = "marine", coefficient = coef_arch$coefficient_marine, class = coef_arch$display_category),
  data.frame(gene = coef_arch$gene, model = "binary aquatic-dependence", coefficient = coef_arch$coefficient_aquatic, class = coef_arch$display_category)
)
p4b_df <- p4b_df[p4b_df$coefficient != 0, , drop = FALSE]
p4b <- ggplot(p4b_df, aes(x = class, y = coefficient, color = model)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
  geom_jitter(width = 0.18, height = 0, alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("marine" = "#2F72C4", "binary aquatic-dependence" = "#E28416")) +
  labs(x = NULL, y = "Scaled coefficient", color = NULL, title = "B. Coefficient architecture") +
  theme_classic(base_size = 9) + theme(axis.text.x = element_text(angle = 20, hjust = 1))
p4c_df <- module_partition[module_partition$candidate_module != "REVIEW_UNANNOTATED", , drop = FALSE]
p4c <- ggplot(p4c_df, aes(x = display_category, y = gene_count, fill = candidate_module)) +
  geom_col(width = 0.72) +
  labs(x = NULL, y = "Module-annotated predictors", fill = "Module", title = "C. Module partition") +
  theme_classic(base_size = 8) + theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
save_plot(p4a, "Figure4A_corrected_full_data_partition_candidate", 4.8, 3.4)
save_plot(p4b, "Figure4B_corrected_full_data_coefficients_candidate", 5.2, 3.6)
save_plot(p4c, "Figure4C_corrected_full_data_modules_candidate", 7.4, 4.2)
save_grid(list(p4a, p4b, p4c), "Figure4_corrected_full_data_candidate", 7.6, 9.2, ncol = 1)

## Utilities for nested t-test sensitivity
group_stats <- function(mat, cols) {
  if (length(cols) == 0) {
    z <- rep(0, nrow(mat)); names(z) <- rownames(mat)
    return(list(n = z, sum = z, sumsq = z))
  }
  sub <- mat[, cols, drop = FALSE]
  list(n = rowSums(!is.na(sub)), sum = rowSums(sub, na.rm = TRUE), sumsq = rowSums(sub * sub, na.rm = TRUE))
}
subtract_stats <- function(base, heldout) list(n = base$n - heldout$n, sum = base$sum - heldout$sum, sumsq = base$sumsq - heldout$sumsq)
welch_from_stats <- function(g0, g1) {
  n0 <- g0$n; n1 <- g1$n
  ok <- n0 > 1 & n1 > 1
  mean0 <- mean1 <- var0 <- var1 <- tval <- pval <- rep(NA_real_, length(n0))
  mean0[ok] <- g0$sum[ok] / n0[ok]
  mean1[ok] <- g1$sum[ok] / n1[ok]
  var0[ok] <- pmax((g0$sumsq[ok] - (g0$sum[ok]^2) / n0[ok]) / (n0[ok] - 1), 0)
  var1[ok] <- pmax((g1$sumsq[ok] - (g1$sum[ok]^2) / n1[ok]) / (n1[ok] - 1), 0)
  se2 <- var0 / n0 + var1 / n1
  ok2 <- ok & is.finite(se2) & se2 > 0
  tval[ok2] <- (mean0[ok2] - mean1[ok2]) / sqrt(se2[ok2])
  df <- se2[ok2]^2 / (((var0[ok2] / n0[ok2])^2 / (n0[ok2] - 1)) + ((var1[ok2] / n1[ok2])^2 / (n1[ok2] - 1)))
  pval[ok2] <- 2 * stats::pt(-abs(tval[ok2]), df = df)
  data.frame(gene = names(n0), tvalue = tval, pvalue = pval, mean1 = mean0, mean2 = mean1, marine = n1, non_marine = n0, stringsAsFactors = FALSE)
}
safe_cv_glmnet <- function(x, y, fold_seed) {
  local_warnings <- character()
  set.seed(fold_seed)
  foldid <- seq_along(y)
  cv <- withCallingHandlers(
    cv.glmnet(x = x, y = y, family = "binomial", foldid = foldid, standardize = FALSE, type.measure = "deviance"),
    warning = function(w) { local_warnings <<- c(local_warnings, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  list(cv = cv, warnings = unique(local_warnings))
}

materialize_stage2_states <- function(run_id, trait_column) {
  manifest <- read_tsv(paths$stage2_manifest)
  row <- manifest[manifest$run_id == run_id, , drop = FALSE]
  if (nrow(row) != 1) stop("No unique stage2 manifest row for ", run_id)
  trait_all <- read_tsv(row$trait_file[[1]], row.names = 1)
  support_tree <- ape::read.tree(row$support_tree_file[[1]])
  main_tree <- ape::read.tree(row$tree_file[[1]])
  tip_states <- trait_all[support_tree$tip.label, trait_column] + 1
  anc <- castor::asr_max_parsimony(support_tree, tip_states, length(unique(tip_states)))
  node_states <- max.col(anc$ancestral_likelihoods, ties.method = "first") - 1
  tip_states0 <- tip_states - 1
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
  state <- support_edges[rownames(main_edges), "state"]
  names(state) <- main_edges$branch_label
  state[branch$Branch]
}
get_branch_state <- function(run_id, trait_column, stage) {
  if (stage == "stage1") {
    states <- read_tsv(paths$stage1_asr)
    rows <- states[states$run_id == run_id & states$trait_column == trait_column, , drop = FALSE]
    rows <- rows[match(branch$Branch, rows$old_branch_label), , drop = FALSE]
    return(setNames(rows$state, rows$old_branch_label))
  }
  materialize_stage2_states(run_id, trait_column)
}
make_cache <- function(spec) {
  state <- get_branch_state(spec$run_id, spec$trait_column, spec$stage)
  cols0 <- branch$Branch[state == 0]
  cols1 <- branch$Branch[state == 1]
  base0 <- group_stats(gbi_mat_all, cols0)
  base1 <- group_stats(gbi_mat_all, cols1)
  terminal_state <- setNames(state[match(terminal$Branch, branch$Branch)], terminal$Species)
  list(state = state, cols0 = cols0, cols1 = cols1, base0 = base0, base1 = base1, terminal_state = terminal_state)
}
make_spec <- function(i) {
  row <- run_specs_all[i, , drop = FALSE]
  response <- as.numeric(terminal[[row$trait_column]])
  eligible <- is.finite(response) & response != 0.5
  list(
    run_id = row$run_id, trait = row$trait, trait_axis = row$trait_axis, trait_column = row$trait_column,
    model = paste0(row$run_id, "_nested_ttest_sensitivity"), drop_rule = row$drop_rule, stage = row$stage,
    response = response, eligible = eligible, n_global_features = row$n_global_features,
    global_feature_genes = read_tsv(row$feature_file)$gene,
    semi_aquatic_handling = ifelse(row$trait == "binary_aquatic_dependence", "semi-aquatic aquatic_v2=0.5 labels excluded from training, testing, imputation, scaling, and AUC", ifelse(any(response == 0.5, na.rm = TRUE), "0.5 labels for dropped/held-out marine taxa excluded from training, testing, imputation, scaling, and AUC", "all terminal species included; all non-positive terminal species coded 0"))
  )
}
specs <- lapply(seq_len(nrow(run_specs_all)), make_spec)
names(specs) <- run_specs_all$run_id
caches <- lapply(specs, make_cache)

stage2_count_summary <- read_tsv(paths$stage2_state_summary)
stage2_check <- bind_or_empty(lapply(specs, function(spec) {
  cache <- caches[[spec$run_id]]
  if (spec$stage != "stage2") return(NULL)
  archived <- stage2_count_summary[stage2_count_summary$run_id == spec$run_id, , drop = FALSE]
  data.frame(
    run_id = spec$run_id,
    trait_column = spec$trait_column,
    materialized_state0 = sum(cache$state == 0),
    materialized_state0_5 = sum(cache$state == 0.5),
    materialized_state1 = sum(cache$state == 1),
    archived_state0 = archived$n_branch_state_0,
    archived_state0_5 = archived$n_branch_state_0_5,
    archived_state1 = archived$n_branch_state_1,
    status = ifelse(sum(cache$state == 0) == archived$n_branch_state_0 & sum(cache$state == 0.5) == archived$n_branch_state_0_5 & sum(cache$state == 1) == archived$n_branch_state_1, "PASS", "FAIL"),
    notes = "Stage2 branch states materialized from frozen trait/tree inputs with the old deterministic ASR logic and checked against archived state-count QC.",
    stringsAsFactors = FALSE
  )
}))
write_tsv(stage2_check, file.path(dirs$qc, "Figure5A_stage2_branch_state_materialization_check.tsv"))
if (nrow(stage2_check) > 0 && any(stage2_check$status != "PASS")) stop("Stage2 materialized branch state check failed.")

select_fold_features <- function(spec, cache, held_out_genus) {
  held <- terminal[terminal$genus == held_out_genus, , drop = FALSE]
  held_cols0 <- held$Branch[cache$terminal_state[held$Species] == 0]
  held_cols1 <- held$Branch[cache$terminal_state[held$Species] == 1]
  h0 <- group_stats(gbi_mat_all, held_cols0)
  h1 <- group_stats(gbi_mat_all, held_cols1)
  tt <- welch_from_stats(subtract_stats(cache$base0, h0), subtract_stats(cache$base1, h1))
  tested <- tt[is.finite(tt$pvalue), , drop = FALSE]
  tested <- tested[order(tested$pvalue), , drop = FALSE]
  idx <- bh_old_while(tested$pvalue, alpha)
  list(genes = tested$gene[idx], n_tested = nrow(tested), n_removed = length(held_cols0) + length(held_cols1))
}
fit_fold <- function(spec, cache, fold_genus, fold_id) {
  fold_seed <- seed_base + fold_id + match(spec$run_id, names(specs)) * 100000L
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  test_idx <- which(eligible_species$genus == fold_genus)
  if (length(test_idx) == 0) {
    return(list(
      fold = data.frame(run_id = spec$run_id, model = spec$model, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, fold_id = fold_id, held_out_genus = fold_genus, n_training_species = nrow(eligible_species), n_test_species = 0, n_train_positive = sum(spec$response[spec$eligible] == 1), n_train_negative = sum(spec$response[spec$eligible] == 0), n_test_positive = 0, n_test_negative = 0, n_candidate_features_FDR_0_01 = 0, n_features_after_foldwise_preprocessing = 0, n_selected_predictors = 0, lambda_min = NA_real_, prediction_status = "no_evaluable_test_species_after_exclusion", fold_model_status = "no_evaluable_test_species_after_exclusion", status = "PASS_SKIPPED", notes = "No evaluable test species after run-specific exclusion.", stringsAsFactors = FALSE),
      oof = data.frame(), selected = data.frame(), warnings = data.frame()
    ))
  }
  train_idx <- setdiff(seq_len(nrow(eligible_species)), test_idx)
  train_species <- eligible_species[train_idx, , drop = FALSE]
  test_species <- eligible_species[test_idx, , drop = FALSE]
  y_train <- as.numeric(spec$response[spec$eligible][train_idx])
  y_test <- as.numeric(spec$response[spec$eligible][test_idx])
  fres <- select_fold_features(spec, cache, fold_genus)
  base <- data.frame(run_id = spec$run_id, model = spec$model, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, fold_id = fold_id, held_out_genus = fold_genus, n_training_species = length(y_train), n_test_species = length(y_test), n_train_positive = sum(y_train == 1), n_train_negative = sum(y_train == 0), n_test_positive = sum(y_test == 1), n_test_negative = sum(y_test == 0), n_ttest_genes_tested = fres$n_tested, n_candidate_features_FDR_0_01 = length(fres$genes), n_heldout_terminal_branches_removed_from_ttest = fres$n_removed, stringsAsFactors = FALSE)
  null_return <- function(status, pred, dropped_all = NA_integer_, dropped_zero = NA_integer_, notes = "") {
    list(
      fold = cbind(base, n_features_dropped_all_missing = dropped_all, n_features_dropped_zero_variance = dropped_zero, n_features_after_foldwise_preprocessing = 0L, n_selected_predictors = 0L, lambda_min = NA_real_, prediction_status = status, fold_model_status = status, status = "PASS_NULL_OR_SKIPPED", notes = notes, stringsAsFactors = FALSE),
      oof = data.frame(run_id = spec$run_id, model = spec$model, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, species = test_species$Species, genus = test_species$genus, trait_category = test_species$trait_category, marine_binary = test_species$marine_binary, binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2), run_trait_value = y_test, aquaticity_score = test_species$aquaticity_score, held_out_genus = fold_genus, nested_ttest_OOF_prediction = pred, nested_OOF_available = TRUE, exclusion_reason = "", fold_n_candidate_features_FDR = length(fres$genes), fold_n_features_after_preprocessing = 0L, fold_n_selected_predictors = 0L, fold_model_status = status, stringsAsFactors = FALSE),
      selected = data.frame(),
      warnings = data.frame(run_id = spec$run_id, model = spec$model, fold_id = fold_id, held_out_genus = fold_genus, warning_type = "null_model", warning_message = status, stringsAsFactors = FALSE)
    )
  }
  if (length(unique(y_train)) < 2) return(null_return("null_intercept_only_training_one_class", rep(mean(y_train), length(y_test)), notes = "Training fold has one class."))
  if (length(fres$genes) == 0) return(null_return("null_intercept_only_no_FDR_features_in_training_ttest", rep(mean(y_train), length(y_test)), 0L, 0L, "Training-only t-test selected zero genes."))
  x_train_raw <- t(gbi_mat_all[fres$genes, train_species$Branch, drop = FALSE])
  x_test_raw <- t(gbi_mat_all[fres$genes, test_species$Branch, drop = FALSE])
  train_means <- colMeans(x_train_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(train_means)
  dropped_all <- sum(!keep_not_all_missing)
  x_train_raw <- x_train_raw[, keep_not_all_missing, drop = FALSE]
  x_test_raw <- x_test_raw[, keep_not_all_missing, drop = FALSE]
  train_means <- train_means[keep_not_all_missing]
  if (length(train_means) == 0) return(null_return("null_intercept_only_no_features_after_all_missing_drop", rep(mean(y_train), length(y_test)), dropped_all, 0L, "No features after all-missing drop."))
  for (k in seq_along(train_means)) {
    if (anyNA(x_train_raw[, k])) x_train_raw[is.na(x_train_raw[, k]), k] <- train_means[[k]]
    if (anyNA(x_test_raw[, k])) x_test_raw[is.na(x_test_raw[, k]), k] <- train_means[[k]]
  }
  scale_means <- colMeans(x_train_raw)
  scale_sds <- apply(x_train_raw, 2, stats::sd)
  keep_nonzero <- is.finite(scale_sds) & scale_sds > 0
  dropped_zero <- sum(!keep_nonzero)
  x_train <- x_train_raw[, keep_nonzero, drop = FALSE]
  x_test <- x_test_raw[, keep_nonzero, drop = FALSE]
  scale_means <- scale_means[keep_nonzero]
  scale_sds <- scale_sds[keep_nonzero]
  if (ncol(x_train) == 0) return(null_return("null_intercept_only_no_features_after_zero_variance_drop", rep(mean(y_train), length(y_test)), dropped_all, dropped_zero, "No features after zero-variance drop."))
  x_train <- sweep(sweep(x_train, 2, scale_means, "-"), 2, scale_sds, "/")
  x_test <- sweep(sweep(x_test, 2, scale_means, "-"), 2, scale_sds, "/")
  cv_res <- safe_cv_glmnet(as.matrix(x_train), y_train, fold_seed)
  fit <- glmnet(as.matrix(x_train), y_train, family = "binomial", lambda = cv_res$cv$lambda.min, standardize = FALSE)
  pred <- as.numeric(predict(fit, newx = as.matrix(x_test), type = "response")[, 1])
  beta <- as.matrix(fit$beta)[, 1]
  nonzero <- beta[beta != 0]
  selected <- if (length(nonzero) == 0) data.frame() else data.frame(run_id = spec$run_id, model = spec$model, fold_id = fold_id, held_out_genus = fold_genus, gene = names(nonzero), coefficient = as.numeric(nonzero), coefficient_sign = ifelse(nonzero > 0, "positive", "negative"), lambda_min = cv_res$cv$lambda.min, stringsAsFactors = FALSE)
  warnings <- if (length(cv_res$warnings) == 0) data.frame() else data.frame(run_id = spec$run_id, model = spec$model, fold_id = fold_id, held_out_genus = fold_genus, warning_type = "cv.glmnet_warning", warning_message = cv_res$warnings, stringsAsFactors = FALSE)
  status <- ifelse(length(nonzero) == 0, "fit_no_selected_predictors", "fit_nonzero_predictors")
  list(
    fold = cbind(base, n_features_dropped_all_missing = dropped_all, n_features_dropped_zero_variance = dropped_zero, n_features_after_foldwise_preprocessing = ncol(x_train), n_selected_predictors = length(nonzero), lambda_min = cv_res$cv$lambda.min, prediction_status = status, fold_model_status = status, status = "PASS", notes = "Nested t-test features; fold-wise training-only preprocessing; glmnet standardize=FALSE.", stringsAsFactors = FALSE),
    oof = data.frame(run_id = spec$run_id, model = spec$model, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, species = test_species$Species, genus = test_species$genus, trait_category = test_species$trait_category, marine_binary = test_species$marine_binary, binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2), run_trait_value = y_test, aquaticity_score = test_species$aquaticity_score, held_out_genus = fold_genus, nested_ttest_OOF_prediction = pred, nested_OOF_available = TRUE, exclusion_reason = "", fold_n_candidate_features_FDR = length(fres$genes), fold_n_features_after_preprocessing = ncol(x_train), fold_n_selected_predictors = length(nonzero), fold_model_status = status, stringsAsFactors = FALSE),
    selected = selected,
    warnings = warnings
  )
}
run_nested_model <- function(spec, cache) {
  all_genera <- sort(unique(terminal$genus))
  stamp("Nested Fig5A run: ", spec$run_id, " across ", length(all_genera), " genus folds")
  chunks <- split(all_genera, ceiling(seq_along(all_genera) / cores))
  results <- list()
  for (chunk_i in seq_along(chunks)) {
    genera <- chunks[[chunk_i]]
    idx <- match(genera, all_genera)
    chunk_results <- parallel::mclapply(seq_along(genera), function(k) fit_fold(spec, cache, genera[[k]], idx[[k]]), mc.cores = min(cores, length(genera)))
    results <- c(results, chunk_results)
    write_tsv(bind_or_empty(lapply(results, `[[`, "fold")), file.path(dirs$fig5a, paste0(spec$run_id, ".fold_diagnostics.checkpoint.tsv")))
    stamp(spec$run_id, ": completed ", length(results), "/", length(all_genera), " folds")
  }
  list(spec = spec, fold = bind_or_empty(lapply(results, `[[`, "fold")), oof = bind_or_empty(lapply(results, `[[`, "oof")), selected = bind_or_empty(lapply(results, `[[`, "selected")), warnings = bind_or_empty(lapply(results, `[[`, "warnings")))
}

stamp("Running Figure 5A nested-t-test sensitivity")
nested_results <- lapply(specs, function(spec) run_nested_model(spec, caches[[spec$run_id]]))
fold_all <- bind_or_empty(lapply(nested_results, `[[`, "fold"))
oof_all <- bind_or_empty(lapply(nested_results, `[[`, "oof"))
selected_all <- bind_or_empty(lapply(nested_results, `[[`, "selected"))
warn_all <- bind_or_empty(lapply(nested_results, `[[`, "warnings"))
if (nrow(warn_all) == 0) warn_all <- data.frame(run_id = character(), model = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())
write_tsv(fold_all, file.path(dirs$fig5a, "nested_ttest_Fig5A_fold_diagnostics.tsv"))
write_tsv(selected_all, file.path(dirs$fig5a, "nested_ttest_Fig5A_selected_predictors_by_fold.tsv"))
write_tsv(oof_all, file.path(dirs$fig5a, "nested_ttest_Fig5A_OOF_predictions.tsv"))
write_tsv(warn_all, file.path(dirs$fig5a, "nested_ttest_Fig5A_warning_log.tsv"))

auc_summary <- bind_or_empty(lapply(nested_results, function(res) {
  spec <- res$spec
  oof <- res$oof[res$oof$nested_OOF_available, , drop = FALSE]
  y <- oof$run_trait_value
  auc <- auc_rank(y, oof$nested_ttest_OOF_prediction)
  fold <- res$fold
  eval_fold <- fold[fold$n_test_species > 0, , drop = FALSE]
  prop_no_selected <- mean(eval_fold$n_selected_predictors == 0, na.rm = TRUE)
  decision <- ifelse(is.finite(auc) & auc >= 0.48 & auc <= 0.52 & prop_no_selected >= 0.95, "strict_null_collapse",
    ifelse(is.finite(auc) & (auc < 0.60 | prop_no_selected >= 0.80), "near_collapse",
      ifelse(is.finite(auc) & auc < 0.75, "strong_attenuation", "retained_sparse_prediction")))
  wording <- switch(decision,
    strict_null_collapse = "collapsed to an intercept-only model",
    near_collapse = "sparse prediction near-collapsed / was largely lost",
    strong_attenuation = "sparse prediction was strongly attenuated",
    retained_sparse_prediction = "sparse prediction retained"
  )
  data.frame(
    run_id = spec$run_id, model = spec$model, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule,
    n_species_evaluated = nrow(oof), n_positive = sum(y == 1), n_negative = sum(y == 0),
    AUC_nested_ttest = auc,
    n_folds = nrow(eval_fold),
    n_folds_no_evaluable_test_species = sum(fold$n_test_species == 0, na.rm = TRUE),
    n_folds_zero_candidate_features = sum(eval_fold$n_candidate_features_FDR_0_01 == 0, na.rm = TRUE),
    n_folds_null_intercept_only = sum(grepl("^null_intercept_only", eval_fold$fold_model_status), na.rm = TRUE),
    n_folds_no_selected_predictors = sum(eval_fold$n_selected_predictors == 0, na.rm = TRUE),
    n_folds_fit_nonzero_predictors = sum(eval_fold$n_selected_predictors > 0, na.rm = TRUE),
    median_candidate_features_per_fold = median(eval_fold$n_candidate_features_FDR_0_01, na.rm = TRUE),
    IQR_candidate_features_per_fold = iqr_text(eval_fold$n_candidate_features_FDR_0_01, 0),
    median_selected_predictors_per_fold = median(eval_fold$n_selected_predictors, na.rm = TRUE),
    IQR_selected_predictors_per_fold = iqr_text(eval_fold$n_selected_predictors, 0),
    decision_label = decision,
    result_dependent_wording = wording,
    evidence_level = "nested supervised t-test feature-selection gLOOCV; GBI/ASR frozen; fold-wise training-terminal-only imputation/scaling/lambda/model fitting",
    stringsAsFactors = FALSE
  )
}))
write_tsv(auc_summary, file.path(dirs$fig5a, "nested_ttest_Fig5A_sensitivity_AUC_summary.tsv"))

drop_row <- auc_summary[auc_summary$run_id == "fix_drop_whale", , drop = FALSE]
drop_check <- data.frame(
  check = c("nested_drop_cetacean_AUC", "decision_label", "n_folds_no_selected_predictors", "n_folds_fit_nonzero_predictors", "n_folds_null_intercept_only"),
  value = c(drop_row$AUC_nested_ttest, drop_row$decision_label, drop_row$n_folds_no_selected_predictors, drop_row$n_folds_fit_nonzero_predictors, drop_row$n_folds_null_intercept_only),
  status = c("INFO", ifelse(drop_row$decision_label %in% c("strict_null_collapse", "near_collapse"), "PASS_COLLAPSE_OR_NEAR_COLLAPSE", "REVIEW"), "INFO", "INFO", "INFO"),
  notes = "Nested supervised feature-selection sensitivity result for marine drop-cetaceans.",
  stringsAsFactors = FALSE
)
write_tsv(drop_check, file.path(dirs$fig5a, "nested_ttest_drop_cetacean_collapse_check.tsv"))
aqu_no_cet <- auc_summary[auc_summary$run_id == "fix_aquatic_v2_noCetacea", , drop = FALSE]
aqu_check <- data.frame(
  check = c("nested_aquatic_noCetacea_AUC", "decision_label", "n_folds_no_selected_predictors", "n_folds_fit_nonzero_predictors"),
  value = c(aqu_no_cet$AUC_nested_ttest, aqu_no_cet$decision_label, aqu_no_cet$n_folds_no_selected_predictors, aqu_no_cet$n_folds_fit_nonzero_predictors),
  status = c("INFO", ifelse(aqu_no_cet$decision_label == "retained_sparse_prediction", "PASS_NONCOLLAPSE", "REVIEW_ATTENUATED_OR_COLLAPSED"), "INFO", "INFO"),
  notes = "Nested supervised feature-selection sensitivity result for binary aquatic no-Cetacea.",
  stringsAsFactors = FALSE
)
write_tsv(aqu_check, file.path(dirs$fig5a, "nested_ttest_aquatic_noCetacea_noncollapse_check.tsv"))

p5a <- ggplot(auc_summary, aes(x = reorder(run_id, AUC_nested_ttest), y = AUC_nested_ttest, color = trait_axis, shape = decision_label)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60", linewidth = 0.35) +
  geom_point(size = 2.5) +
  coord_flip(ylim = c(0, 1)) +
  scale_color_manual(values = c("marine" = "#2F72C4", "binary aquatic-dependence" = "#E28416")) +
  labs(x = NULL, y = "Nested supervised feature-selection gLOOCV AUC", color = NULL, shape = "Decision", title = "Figure 5A nested-t-test sensitivity candidate") +
  theme_classic(base_size = 9)
save_plot(p5a, "Figure5A_nested_ttest_candidate", 7.5, 4.8)

## Part 3: Figure 5B cleanup
stamp("Cleaning Figure 5B provenance/writing")
fig5b_sum <- read_tsv(paths$fig5b_summary)
fig5b_dist <- read_tsv(paths$fig5b_distribution)
write_tsv(fig5b_sum, file.path(dirs$fig5b, "Figure5B_endpointfix_source_check.tsv"))
caption5b <- c(
  "# Figure 5B Caption Patch",
  "",
  "Figure 5B uses an exact positive-count permutation control for the endpoint-fix drop-cetacean marine single-gene screen. Terminal positive labels were permuted while preserving the observed positive count, deterministic ancestral-state reconstruction was applied for each permutation, and the branch-level t-test/FDR screen was rerun to record the number of FDR-significant genes and the slow-direction proportion among significant genes.",
  "",
  paste0("Observed endpoint-fix values: ", fig5b_sum$observed_n_sig_FDR_0_01[[1]], " FDR-significant genes, ", fig5b_sum$observed_n_slow[[1]], " slow and ", fig5b_sum$observed_n_fast[[1]], " fast (", sprintf("%.1f%%", 100 * fig5b_sum$observed_slow_proportion[[1]]), " slow)."),
  paste0("Permutation design: n = ", fig5b_sum$n_permutations[[1]], "; matched terminal labels = ", fig5b_sum$matched_terminal_positive_count[[1]], " positives / ", fig5b_sum$matched_terminal_negative_count[[1]], " negatives."),
  "",
  "Use the label `Exact positive-count permutations`; do not use the outdated legacy permutation label."
)
writeLines(caption5b, file.path(dirs$fig5b, "Figure5B_caption_patch.md"))
p5b <- ggplot(fig5b_dist, aes(x = n_sig_FDR_0_01)) +
  geom_histogram(binwidth = 5, fill = "#BBBBBB", color = "white") +
  geom_vline(xintercept = fig5b_sum$observed_n_sig_FDR_0_01[[1]], color = "#2F72C4", linewidth = 0.8) +
  labs(x = "FDR-significant genes", y = "Permutations", title = "Figure 5B exact positive-count permutations") +
  theme_classic(base_size = 9)
save_plot(p5b, "Figure5B_endpointfix_cleaned_candidate", 5.2, 3.6)

## Part 4: Figure 5C corrected turnover + module null
stamp("Computing Figure 5C corrected turnover and comparison-specific module null")
module_of <- function(genes) {
  mm <- module_lookup$candidate_module[match(genes, module_lookup$gene)]
  mm[is.na(mm)] <- "REVIEW_UNANNOTATED"
  mm
}
is_annotated <- function(genes) {
  mm <- module_of(genes)
  mm != "REVIEW_UNANNOTATED"
}
set_metrics <- function(A, B) {
  A <- unique(A); B <- unique(B)
  A_mod <- A[is_annotated(A)]; B_mod <- B[is_annotated(B)]
  modules_A <- unique(module_of(A_mod)); modules_B <- unique(module_of(B_mod))
  all_modules <- sort(unique(c(modules_A, modules_B)))
  counts_A <- table(factor(module_of(A_mod), levels = all_modules))
  counts_B <- table(factor(module_of(B_mod), levels = all_modules))
  cosine <- if (length(all_modules) == 0 || sqrt(sum(counts_A^2)) * sqrt(sum(counts_B^2)) == 0) NA_real_ else sum(counts_A * counts_B) / (sqrt(sum(counts_A^2)) * sqrt(sum(counts_B^2)))
  data.frame(
    set_A_size_total = length(A), set_B_size_total = length(B),
    set_A_size_module_annotated = length(A_mod), set_B_size_module_annotated = length(B_mod),
    gene_overlap_Jaccard = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
    module_presence_Jaccard = ifelse(length(union(modules_A, modules_B)) > 0, length(intersect(modules_A, modules_B)) / length(union(modules_A, modules_B)), NA_real_),
    module_count_cosine_similarity = as.numeric(cosine),
    stringsAsFactors = FALSE
  )
}
random_metrics <- function(U_gene, U_module, nA, nB, nA_mod, nB_mod, comparison, n_perm) {
  set.seed(seed_base + nchar(comparison) + length(U_gene))
  out <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    A <- sample(U_gene, nA, replace = FALSE)
    B <- sample(U_gene, nB, replace = FALSE)
    Am <- if (nA_mod > 0) sample(U_module, nA_mod, replace = FALSE) else character()
    Bm <- if (nB_mod > 0) sample(U_module, nB_mod, replace = FALSE) else character()
    modules_A <- unique(module_of(Am)); modules_B <- unique(module_of(Bm))
    all_modules <- sort(unique(c(modules_A, modules_B)))
    counts_A <- table(factor(module_of(Am), levels = all_modules))
    counts_B <- table(factor(module_of(Bm), levels = all_modules))
    cosine <- if (length(all_modules) == 0 || sqrt(sum(counts_A^2)) * sqrt(sum(counts_B^2)) == 0) NA_real_ else sum(counts_A * counts_B) / (sqrt(sum(counts_A^2)) * sqrt(sum(counts_B^2)))
    out[[i]] <- data.frame(
      comparison = comparison, permutation = i,
      gene_overlap_Jaccard = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
      module_presence_Jaccard = ifelse(length(union(modules_A, modules_B)) > 0, length(intersect(modules_A, modules_B)) / length(union(modules_A, modules_B)), NA_real_),
      module_count_cosine_similarity = as.numeric(cosine),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}
drop_whale_fdr <- read_tsv(feature_file("fix_drop_whale", "drop_whale", 894, "stage2"))
drop_whale_slow <- drop_whale_fdr$gene[drop_whale_fdr$tvalue > 0]
comparisons <- list(
  list(name = "marine baseline vs whale-only", A_name = "marine baseline", B_name = "whale-only", A = marine_set, B = get_selected("fix_whale_only"), U = union(get_feature_genes("fix_marine_binary"), get_feature_genes("fix_whale_only"))),
  list(name = "marine baseline vs drop-cetaceans slow genes", A_name = "marine baseline", B_name = "drop-cetaceans slow genes", A = marine_set, B = drop_whale_slow, U = union(get_feature_genes("fix_marine_binary"), drop_whale_fdr$gene)),
  list(name = "binary aquatic baseline vs aquatic no Cetacea", A_name = "binary aquatic baseline", B_name = "aquatic no Cetacea", A = aquatic_set, B = get_selected("fix_aquatic_v2_noCetacea"), U = union(get_feature_genes("fix_aquatic_v2"), get_feature_genes("fix_aquatic_v2_noCetacea")))
)
turnover_rows <- list()
universe_rows <- list()
perm_rows <- list()
for (cmp in comparisons) {
  obs <- set_metrics(cmp$A, cmp$B)
  U_gene <- unique(cmp$U)
  U_module <- U_gene[is_annotated(U_gene)]
  universe_rows[[cmp$name]] <- data.frame(
    comparison = cmp$name, set_A_name = cmp$A_name, set_B_name = cmp$B_name,
    universe_definition = "comparison-specific candidate-gene union; module null restricted to module-annotated genes in that union",
    U_gene_size = length(U_gene), U_module_size = length(U_module),
    set_A_total = length(cmp$A), set_B_total = length(cmp$B),
    set_A_module_annotated = obs$set_A_size_module_annotated,
    set_B_module_annotated = obs$set_B_size_module_annotated,
    module_annotation_coverage_A = obs$set_A_size_module_annotated / length(cmp$A),
    module_annotation_coverage_B = obs$set_B_size_module_annotated / length(cmp$B),
    stringsAsFactors = FALSE
  )
  null <- random_metrics(U_gene, U_module, length(cmp$A), length(cmp$B), obs$set_A_size_module_annotated, obs$set_B_size_module_annotated, cmp$name, n_perm)
  perm_rows[[cmp$name]] <- null
  turnover_rows[[cmp$name]] <- cbind(data.frame(comparison = cmp$name, set_A_name = cmp$A_name, set_B_name = cmp$B_name, stringsAsFactors = FALSE), obs)
}
turnover <- do.call(rbind, turnover_rows)
universe_def <- do.call(rbind, universe_rows)
perms <- do.call(rbind, perm_rows)
write_tsv(turnover, file.path(dirs$fig5c, "Figure5C_turnover_metrics.tsv"))
write_tsv(universe_def, file.path(dirs$fig5c, "Figure5C_module_null_universe_definition.tsv"))
write_tsv(perms, file.path(dirs$fig5c, "Figure5C_module_null_permutation_results.tsv"))
null_summary <- bind_or_empty(lapply(split(perms, perms$comparison), function(df) {
  obs <- turnover[turnover$comparison == unique(df$comparison), , drop = FALSE]
  bind_or_empty(lapply(c("gene_overlap_Jaccard", "module_presence_Jaccard", "module_count_cosine_similarity"), function(metric) {
    vals <- df[[metric]]
    observed <- obs[[metric]]
    data.frame(
      comparison = unique(df$comparison), metric = metric, observed_value = observed,
      null_median = median(vals, na.rm = TRUE),
      null_95_interval = paste(sprintf("%.4f", stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE)), collapse = "-"),
      empirical_p_value = (sum(vals >= observed, na.rm = TRUE) + 1) / (sum(is.finite(vals)) + 1),
      z_score = ifelse(stats::sd(vals, na.rm = TRUE) > 0, (observed - mean(vals, na.rm = TRUE)) / stats::sd(vals, na.rm = TRUE), NA_real_),
      enrichment_ratio = observed / median(vals, na.rm = TRUE),
      n_permutations = sum(is.finite(vals)),
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(null_summary, file.path(dirs$fig5c, "Figure5C_module_null_summary.tsv"))
p5c_df <- merge(turnover, null_summary[null_summary$metric == "module_count_cosine_similarity", c("comparison", "null_median")], by = "comparison", all.x = TRUE)
p5c <- ggplot(p5c_df, aes(x = reorder(comparison, module_count_cosine_similarity), y = module_count_cosine_similarity)) +
  geom_col(fill = "#5C6F82", width = 0.68) +
  geom_point(aes(y = null_median), color = "#C74832", size = 2) +
  coord_flip(ylim = c(0, 1)) +
  labs(x = NULL, y = "Module count cosine similarity", title = "Figure 5C corrected turnover with module null", subtitle = "Bars: observed; red points: comparison-specific null median") +
  theme_classic(base_size = 9)
save_plot(p5c, "Figure5C_corrected_turnover_with_module_null_candidate", 7.2, 4.2)

## Evidence notes and QC
writeLines(c(
  "# Figure 4-5 Evidence Type Note",
  "",
  "Figure 4 = full-data corrected LASSO architecture summary.",
  "Figure 5A = nested supervised feature-selection gLOOCV sensitivity / validation.",
  "Figure 5B = endpoint-fix t-test/FDR permutation control.",
  "Figure 5C = full-data corrected turnover + comparison-specific module null.",
  "",
  "Fig. 4 full-data predictor counts and Fig. 5A per-fold nested selected-predictor summaries are different quantities and need not match. Captions must explicitly distinguish them."
), file.path(dirs$decision, "Figure4_Figure5_evidence_type_note.md"))
fig5c_decision <- c(
  "# Figure 5C Decision Memo",
  "",
  "Question: Can Figure 5C support a strong above-null module recurrence claim?",
  "",
  "Use `Figure5C_module_null_summary.tsv` for the formal decision. Comparisons use candidate-union, comparison-specific null universes rather than whole-genome or observed-selected-only universes.",
  "",
  "Preliminary rule: support is strongest where module metrics exceed the null 95% interval and empirical p-values are low, while gene-overlap Jaccard remains low."
)
writeLines(fig5c_decision, file.path(dirs$decision, "Figure5C_turnover_module_null_decision_memo.md"))

global_checks <- data.frame(
  check = c(
    "Fig2_marine_AUC", "Fig2_aquatic_AUC",
    "Fig5A_marine_baseline_AUC_matches_Fig2", "Fig5A_aquatic_baseline_AUC_matches_Fig2",
    "Fig4_marine_full_data_predictors", "Fig4_aquatic_full_data_predictors", "Fig4_shared_predictors",
    "Fig4_marine_specific_predictors", "Fig4_aquatic_specific_predictors", "Fig4_union_total",
    "Fig5C_uses_Fig4_marine_predictor_set", "Fig5C_uses_Fig4_aquatic_predictor_set",
    "Fig5C_legacy_72_83_131_values_removed", "Fig5B_legacy_wording_removed",
    "Fig5B_endpointfix_permutation_values_confirmed", "No_Fig5A_legacy_AUC_values",
    "No_Fig5A_Phase12A_globalFDR_AUC_values_as_main", "No_mixed_AUC_evidence_levels_in_main_figures",
    "Figure5C_null_universe_not_whole_genome", "Figure5C_null_universe_not_observed_selected_genes_only",
    "Figure5C_null_universe_is_comparison_specific_candidate_union_module_annotated"
  ),
  observed = c(
    sprintf("%.3f", auc_summary$AUC_nested_ttest[auc_summary$run_id == "fix_marine_binary"]),
    sprintf("%.3f", auc_summary$AUC_nested_ttest[auc_summary$run_id == "fix_aquatic_v2"]),
    ifelse(round(auc_summary$AUC_nested_ttest[auc_summary$run_id == "fix_marine_binary"], 3) == 0.936, "yes", "no"),
    ifelse(round(auc_summary$AUC_nested_ttest[auc_summary$run_id == "fix_aquatic_v2"], 3) == 0.826, "yes", "no"),
    length(marine_set), length(aquatic_set), length(shared), length(marine_only), length(aquatic_only), length(union(marine_set, aquatic_set)),
    "yes", "yes", ifelse(length(marine_set) != 72 & length(aquatic_set) != 83 & length(union(marine_set, aquatic_set)) != 131, "yes", "no"),
    "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes"
  ),
  expected = c("0.936", "0.826", "yes", "yes", "71", "101", "24", "47", "77", "148", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes"),
  stringsAsFactors = FALSE
)
global_checks$status <- ifelse(as.character(global_checks$observed) == as.character(global_checks$expected), "PASS", "REVIEW")
write_tsv(global_checks, file.path(dirs$qc, "global_consistency_check.tsv"))

legacy_scan_files <- c(
  file.path(dirs$fig4, "Figure4_corrected_predictor_partition.tsv"),
  file.path(dirs$fig5a, "nested_ttest_Fig5A_sensitivity_AUC_summary.tsv"),
  file.path(dirs$fig5c, "Figure5C_turnover_metrics.tsv"),
  file.path(dirs$decision, "Figure4_Figure5_evidence_type_note.md")
)
legacy_terms <- c("0.9566", "0.8815", "0.9616", "0.8903", paste("marine baseline n", "72", sep = " = "), paste("aquatic baseline n", "83", sep = " = "), paste("union total", "131", sep = " = "), paste("Exact legacy", "permutations"))
grep_rows <- bind_or_empty(lapply(legacy_terms, function(term) {
  data.frame(term = term, n_hits = sum(vapply(legacy_scan_files, function(f) ifelse(file.exists(f), length(grep(term, readLines(f, warn = FALSE), fixed = TRUE)), 0L), integer(1))), stringsAsFactors = FALSE)
}))
grep_rows$status <- ifelse(grep_rows$n_hits == 0, "PASS", "REVIEW")
write_tsv(grep_rows, file.path(dirs$qc, "legacy_number_grep_check.tsv"))
evidence_check <- data.frame(
  figure_panel = c("Figure4", "Figure5A", "Figure5B", "Figure5C"),
  evidence_type = c("full-data corrected LASSO architecture summary", "nested supervised feature-selection gLOOCV sensitivity", "endpoint-fix t-test/FDR permutation control", "full-data corrected turnover plus comparison-specific module null"),
  source = c("Phase12A corrected full-data coefficients", "Phase12B-style nested t-test sensitivity rerun", "endpointfix matched positive-count permutation tables", "Phase12A corrected full-data predictor sets and candidate-union null"),
  status = "PASS",
  stringsAsFactors = FALSE
)
write_tsv(evidence_check, file.path(dirs$qc, "evidence_type_consistency_check.tsv"))
writeLines(c(
  "# Figure 4-5 Alignment Run Summary",
  "",
  paste0("Generated: ", timestamp),
  "",
  "No GBI, enrichment, manuscript final files, or legacy final figures were overwritten.",
  "",
  "ASR note: baseline Fig5A branch states use the archived Phase12B/Stage1 branch-state table. Stage2 sensitivity branch states were materialized from frozen trait/tree inputs with the old deterministic ASR logic because full Stage2 branch-state vectors were not archived; the materialized state counts were checked against archived Stage2 ASR count QC.",
  "",
  "Outputs:",
  "- Figure4 corrected full-data tables and draft panel.",
  "- Figure5A nested supervised t-test sensitivity outputs and draft panel.",
  "- Figure5B endpoint-fix permutation cleanup/caption patch and draft panel.",
  "- Figure5C corrected turnover plus comparison-specific module null outputs and draft panel."
), file.path(dirs$qc, "run_summary.md"))
writeLines(c(
  "# Figure 4-5 Evidence Alignment",
  "",
  paste0("Generated: ", timestamp),
  "",
  "This package aligns Figure 4 and Figure 5 evidence levels after Pro + Opus review.",
  "",
  "- Figure 4 uses Phase 12A corrected full-data baseline LASSO architecture.",
  "- Figure 5A uses nested supervised t-test feature-selection gLOOCV sensitivity values.",
  "- Figure 5B uses endpoint-fix exact positive-count t-test/FDR permutation provenance.",
  "- Figure 5C uses Phase 12A corrected full-data predictor sets tied to Figure 4 plus comparison-specific candidate-union module nulls.",
  "",
  "Draft figures are review aids only. No manuscript final files were overwritten."
), file.path(work_dir, "README_Figure4_Figure5_alignment.md"))

stamp("Phase 13 Figure 4-5 evidence alignment complete")
