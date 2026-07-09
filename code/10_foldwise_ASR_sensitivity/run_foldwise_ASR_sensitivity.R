#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(castor)
  library(glmnet)
  library(parallel)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "10_reviewer_risk_controls", "15_foldwise_ASR_sensitivity")
code_dir <- file.path(work_dir, "code")
qc_dir <- file.path(work_dir, "qc")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(code_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}
mode <- get_arg("--mode", "full")
cores <- as.integer(get_arg("--cores", min(12L, max(1L, parallel::detectCores() - 1L))))
seed_base <- as.integer(get_arg("--seed", 20260601L))
alpha <- as.numeric(get_arg("--alpha", 0.01))
if (!mode %in% c("smoke", "full")) stop("--mode must be smoke or full")
if (!is.finite(cores) || cores < 1L) cores <- 1L
script_start <- Sys.time()

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  tree = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal302.anno.nwk"),
  support_tree = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal302.anno.BL_support.nwk"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  frozen_asr_branch_states = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "ancestral_state_plots", "endpointfix_stage1_ancestral_state_branch_table.tsv"),
  phase12b_auc = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_AUC_summary.tsv"),
  phase12b_oof = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_OOF_predictions.tsv"),
  phase12b_fold = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_fold_candidate_gene_counts.tsv"),
  phase12b_selected = file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness", "nested_ttest_baseline", "nested_ttest_fold_selected_predictors.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

read_tsv <- function(file, ...) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_tsv <- function(x, file) utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
stamp <- function(...) message("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", ...)

auc_rank <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  n_pos <- sum(y == 1)
  n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

iqr_text <- function(x, digits = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  sprintf(paste0("%.", digits, "f-%.", digits, "f"), qs[[1]], qs[[2]])
}

safe_cv_glmnet <- function(x, y, fold_seed) {
  local_warnings <- character()
  set.seed(fold_seed)
  foldid <- seq_along(y)
  cv <- withCallingHandlers(
    cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      foldid = foldid,
      standardize = FALSE,
      type.measure = "deviance"
    ),
    warning = function(w) {
      local_warnings <<- c(local_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(cv = cv, warnings = unique(local_warnings))
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

group_stats <- function(mat, cols) {
  if (length(cols) == 0) {
    z <- rep(0, nrow(mat))
    names(z) <- rownames(mat)
    return(list(n = z, sum = z, sumsq = z))
  }
  sub <- mat[, cols, drop = FALSE]
  n <- rowSums(!is.na(sub))
  sumv <- rowSums(sub, na.rm = TRUE)
  sumsq <- rowSums(sub * sub, na.rm = TRUE)
  list(n = n, sum = sumv, sumsq = sumsq)
}

welch_from_stats <- function(g0, g1) {
  n0 <- g0$n
  n1 <- g1$n
  ok <- n0 > 1 & n1 > 1
  mean0 <- rep(NA_real_, length(n0))
  mean1 <- rep(NA_real_, length(n0))
  var0 <- rep(NA_real_, length(n0))
  var1 <- rep(NA_real_, length(n0))
  tval <- rep(NA_real_, length(n0))
  pval <- rep(NA_real_, length(n0))
  mean0[ok] <- g0$sum[ok] / n0[ok]
  mean1[ok] <- g1$sum[ok] / n1[ok]
  var0[ok] <- pmax((g0$sumsq[ok] - (g0$sum[ok] * g0$sum[ok]) / n0[ok]) / (n0[ok] - 1), 0)
  var1[ok] <- pmax((g1$sumsq[ok] - (g1$sum[ok] * g1$sum[ok]) / n1[ok]) / (n1[ok] - 1), 0)
  se2 <- var0 / n0 + var1 / n1
  ok2 <- ok & is.finite(se2) & se2 > 0
  tval[ok2] <- (mean0[ok2] - mean1[ok2]) / sqrt(se2[ok2])
  df_num <- se2[ok2] * se2[ok2]
  df_den <- ((var0[ok2] / n0[ok2])^2 / (n0[ok2] - 1)) +
    ((var1[ok2] / n1[ok2])^2 / (n1[ok2] - 1))
  df <- df_num / df_den
  pval[ok2] <- 2 * stats::pt(-abs(tval[ok2]), df = df)
  data.frame(
    gene = names(n0),
    tvalue = tval,
    pvalue = pval,
    mean1 = mean0,
    mean2 = mean1,
    marine = n1,
    non_marine = n0,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

stamp("Reading inputs")
gbi <- read_tsv(paths$gbi)
stopifnot("gene" %in% names(gbi))
rownames(gbi) <- gbi$gene
gbi$gene <- NULL
gbi_mat_all <- as.matrix(gbi)
storage.mode(gbi_mat_all) <- "numeric"
rm(gbi)

branch <- read_tsv(paths$branch)
trait <- read_tsv(paths$trait)
frozen_asr_states <- read_tsv(paths$frozen_asr_branch_states)
phase12b_auc <- read_tsv(paths$phase12b_auc)
phase12b_oof <- read_tsv(paths$phase12b_oof)
phase12b_fold <- read_tsv(paths$phase12b_fold)
phase12b_selected <- read_tsv(paths$phase12b_selected)
support_tree <- ape::read.tree(paths$support_tree)
main_tree <- ape::read.tree(paths$tree)

terminal <- branch[branch$Species != "internal", , drop = FALSE]
terminal$genus <- vapply(strsplit(terminal$Species, "_"), `[`, character(1), 1)
terminal <- merge(terminal, trait, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- terminal[match(branch$Species[branch$Species != "internal"], terminal$Species), , drop = FALSE]
terminal$trait_category <- ifelse(
  terminal$marine_binary == 1, "marine",
  ifelse(terminal$aquatic_v2 == 1, "non-marine aquatic",
    ifelse(terminal$aquatic_v2 == 0.5, "semi-aquatic", "terrestrial")
  )
)
terminal$aquaticity_score <- terminal$aquaticity_score_sum_0_18
if (anyNA(terminal$marine_binary) || anyNA(terminal$aquatic_v2)) {
  stop("Terminal trait mapping contains missing marine_binary or aquatic_v2 values.")
}
if (!all(branch$Branch %in% colnames(gbi_mat_all))) stop("Branch IDs missing from GBI matrix.")
gbi_mat_all <- gbi_mat_all[, branch$Branch, drop = FALSE]

tip_species <- support_tree$tip.label
if (!all(tip_species %in% terminal$Species)) stop("Support tree tips are not all present in terminal branch table.")
trait_by_species <- trait
rownames(trait_by_species) <- trait_by_species$species

main_edges <- data.frame(
  anc = main_tree$edge[, 1],
  offspring = main_tree$edge[, 2],
  branch_label_raw = main_tree$edge.length,
  order = seq_len(nrow(main_tree$edge)),
  stringsAsFactors = FALSE
)
main_edges$old_branch_label <- paste0("B", as.integer(round(main_edges$branch_label_raw)))
main_edges$row_key <- paste(main_edges$anc, main_edges$offspring, sep = "-")

support_edges_template <- data.frame(
  anc = support_tree$edge[, 1],
  offspring = support_tree$edge[, 2],
  support_edge_order = seq_len(nrow(support_tree$edge)),
  stringsAsFactors = FALSE
)
support_edges_template$row_key <- paste(support_edges_template$anc, support_edges_template$offspring, sep = "-")
if (!all(main_edges$row_key %in% support_edges_template$row_key)) {
  stop("Main-tree edges cannot all be matched to support-tree edges by node pair.")
}

make_branch_state_hsp <- function(trait_column, held_out_genus = NULL) {
  tip_state0 <- trait_by_species[tip_species, trait_column]
  tip_input <- tip_state0 + 1
  masked_species <- character()
  if (!is.null(held_out_genus)) {
    masked_species <- terminal$Species[terminal$genus == held_out_genus]
    tip_input[tip_species %in% masked_species] <- NA_real_
  }
  nstates <- length(unique((trait_by_species[, trait_column] + 1)[!is.na(trait_by_species[, trait_column])]))
  hsp <- castor::hsp_max_parsimony(
    tree = support_tree,
    tip_states = tip_input,
    Nstates = nstates,
    transition_costs = "all_equal",
    edge_exponent = 0,
    weight_by_scenarios = TRUE
  )
  if (!isTRUE(hsp$success)) stop("hsp_max_parsimony failed for ", trait_column, " held_out_genus=", held_out_genus)

  all_states <- hsp$states - 1
  # The Stage 1 deterministic ASR table preserves known terminal labels directly.
  known_tip <- !is.na(tip_input)
  all_states[seq_along(tip_state0)][known_tip] <- tip_state0[known_tip]
  names(all_states) <- paste0("B", seq_along(all_states))

  support_edges <- support_edges_template
  support_edges$support_child_node_id <- paste0("B", support_edges$offspring)
  support_edges$state <- all_states[support_edges$support_child_node_id]
  support_edges$node_type <- ifelse(support_edges$offspring <= length(tip_species), "terminal", "internal")
  rownames(support_edges) <- support_edges$row_key
  branch_state <- support_edges[main_edges$row_key, "state"]
  names(branch_state) <- main_edges$old_branch_label
  branch_state <- branch_state[branch$Branch]
  names(branch_state) <- branch$Branch
  list(
    branch_state = branch_state,
    hsp_success = hsp$success,
    n_masked_terminal_labels = length(masked_species),
    masked_species = masked_species
  )
}

model_specs <- list(
  list(
    trait = "marine_binary",
    trait_axis = "marine specialization",
    model = "marine_binary_foldwise_ASR_nested_ttest_baseline",
    frozen_model = "marine_binary_nested_ttest_baseline",
    run_id = "fix_marine_binary",
    trait_column = "marine_binary",
    response_col = "marine_binary",
    eligible = rep(TRUE, nrow(terminal)),
    semi_aquatic_handling = "retained as non-marine/background 0 for marine binary model",
    frozen_auc_public = 0.936
  ),
  list(
    trait = "binary_aquatic_dependence",
    trait_axis = "binary aquatic dependence",
    model = "binary_aquatic_dependence_foldwise_ASR_nested_ttest_baseline",
    frozen_model = "binary_aquatic_dependence_nested_ttest_baseline",
    run_id = "fix_aquatic_v2",
    trait_column = "aquatic_v2",
    response_col = "aquatic_v2",
    eligible = terminal$aquatic_v2 != 0.5,
    semi_aquatic_handling = "semi-aquatic aquatic_v2=0.5 labels excluded from LASSO training/testing/preprocessing/AUC; 0.5 branch states excluded from branch-level binary screen",
    frozen_auc_public = 0.826
  )
)

validate_hsp_against_frozen <- function(spec) {
  h <- make_branch_state_hsp(spec$trait_column, held_out_genus = NULL)
  frozen <- frozen_asr_states[
    frozen_asr_states$run_id == spec$run_id & frozen_asr_states$trait_column == spec$trait_column,
    ,
    drop = FALSE
  ]
  frozen <- frozen[match(branch$Branch, frozen$old_branch_label), , drop = FALSE]
  if (!all(frozen$old_branch_label == branch$Branch)) stop("Frozen ASR order mismatch for ", spec$run_id)
  diff <- frozen$state != h$branch_state
  data.frame(
    trait = spec$trait,
    model = spec$model,
    run_id = spec$run_id,
    hsp_no_missing_reproduces_frozen_ASR = ifelse(sum(diff, na.rm = TRUE) == 0, "yes", "no"),
    n_branch_state_differences_no_missing = sum(diff, na.rm = TRUE),
    n_state0_frozen = sum(frozen$state == 0, na.rm = TRUE),
    n_state0_5_frozen = sum(frozen$state == 0.5, na.rm = TRUE),
    n_state1_frozen = sum(frozen$state == 1, na.rm = TRUE),
    n_state0_hsp = sum(h$branch_state == 0, na.rm = TRUE),
    n_state0_5_hsp = sum(h$branch_state == 0.5, na.rm = TRUE),
    n_state1_hsp = sum(h$branch_state == 1, na.rm = TRUE),
    status = ifelse(sum(diff, na.rm = TRUE) == 0, "PASS", "STOP"),
    notes = "castor::hsp_max_parsimony without masked labels was checked against the existing frozen deterministic ASR branch-state table.",
    stringsAsFactors = FALSE
  )
}

hsp_validation <- do.call(rbind, lapply(model_specs, validate_hsp_against_frozen))
write_tsv(hsp_validation, file.path(qc_dir, "foldwise_ASR_hsp_vs_frozen_validation.tsv"))
if (any(hsp_validation$status != "PASS")) {
  stop("hsp_max_parsimony validation against frozen ASR failed; do not run fold-wise ASR sensitivity.")
}

select_fold_features <- function(spec, branch_state, held_out_genus) {
  held <- terminal[terminal$genus == held_out_genus, , drop = FALSE]
  screen_branches <- names(branch_state)[branch_state %in% c(0, 1)]
  screen_branches <- setdiff(screen_branches, held$Branch)
  cols0 <- screen_branches[branch_state[screen_branches] == 0]
  cols1 <- screen_branches[branch_state[screen_branches] == 1]
  g0 <- group_stats(gbi_mat_all, cols0)
  g1 <- group_stats(gbi_mat_all, cols1)
  tt <- welch_from_stats(g0, g1)
  tt_tested <- tt[is.finite(tt$pvalue), , drop = FALSE]
  tt_tested <- tt_tested[order(tt_tested$pvalue, decreasing = FALSE), , drop = FALSE]
  idx <- bh_old_while(tt_tested$pvalue, alpha = alpha)
  genes <- tt_tested$gene[idx]
  list(
    genes = genes,
    n_ttest_genes_tested = nrow(tt_tested),
    n_branches_screened = length(screen_branches),
    n_screen_state0 = length(cols0),
    n_screen_state1 = length(cols1),
    n_state0_total = sum(branch_state == 0, na.rm = TRUE),
    n_state0_5_total = sum(branch_state == 0.5, na.rm = TRUE),
    n_state1_total = sum(branch_state == 1, na.rm = TRUE),
    n_state_other_total = sum(!(branch_state %in% c(0, 0.5, 1)) & !is.na(branch_state)),
    n_heldout_terminal_branches_removed_from_ttest = sum(held$Branch %in% names(branch_state)[branch_state %in% c(0, 1)])
  )
}

fit_fold <- function(spec, fold_genus, fold_id) {
  fold_seed <- seed_base + fold_id + ifelse(spec$trait == "marine_binary", 100000L, 200000L)
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  test_idx <- which(eligible_species$genus == fold_genus)
  if (length(test_idx) == 0) {
    return(list(
      fold = data.frame(
        trait = spec$trait, model = spec$model, frozen_model = spec$frozen_model,
        fold_id = fold_id, held_out_genus = fold_genus,
        n_training_samples = nrow(eligible_species), n_test_samples = 0,
        n_train_positive = sum(eligible_species[[spec$response_col]] == 1),
        n_train_negative = sum(eligible_species[[spec$response_col]] == 0),
        n_test_positive = 0, n_test_negative = 0,
        n_masked_terminal_labels_for_ASR = length(terminal$Species[terminal$genus == fold_genus]),
        n_ASR_state0_total = NA_integer_, n_ASR_state0_5_total = NA_integer_, n_ASR_state1_total = NA_integer_, n_ASR_state_other_total = NA_integer_,
        n_branches_screened_after_masking_and_heldout_removal = NA_integer_,
        n_screen_state0 = NA_integer_, n_screen_state1 = NA_integer_,
        n_heldout_terminal_branches_removed_from_ttest = NA_integer_,
        n_ttest_genes_tested = NA_integer_, n_fold_candidate_genes = 0L,
        n_features_dropped_all_missing = NA_integer_, n_features_dropped_zero_variance = NA_integer_,
        n_features_used_for_lasso = 0L, n_selected_predictors = 0L,
        lambda_min = NA_real_,
        fold_model_status = "no_evaluable_test_species_after_exclusion",
        fold_seed = fold_seed,
        notes = "No evaluable held-out species for this model after exclusion.",
        stringsAsFactors = FALSE
      ),
      oof = data.frame(),
      selected = data.frame(),
      warnings = data.frame(
        trait = spec$trait, model = spec$model, fold_id = fold_id, held_out_genus = fold_genus,
        warning_type = "fold_skipped",
        warning_message = "no_evaluable_test_species_after_exclusion",
        stringsAsFactors = FALSE
      )
    ))
  }

  train_idx <- setdiff(seq_len(nrow(eligible_species)), test_idx)
  train_species <- eligible_species[train_idx, , drop = FALSE]
  test_species <- eligible_species[test_idx, , drop = FALSE]
  y_train <- as.numeric(train_species[[spec$response_col]])
  y_test <- as.numeric(test_species[[spec$response_col]])
  if (!all(y_train %in% c(0, 1)) || !all(y_test %in% c(0, 1))) stop("Non-binary response in ", spec$model, " fold ", fold_genus)

  asr <- make_branch_state_hsp(spec$trait_column, held_out_genus = fold_genus)
  feature_res <- select_fold_features(spec, asr$branch_state, fold_genus)
  fold_genes <- feature_res$genes
  base_fold <- data.frame(
    trait = spec$trait, model = spec$model, frozen_model = spec$frozen_model,
    fold_id = fold_id, held_out_genus = fold_genus,
    n_training_samples = length(y_train), n_test_samples = length(y_test),
    n_train_positive = sum(y_train == 1), n_train_negative = sum(y_train == 0),
    n_test_positive = sum(y_test == 1), n_test_negative = sum(y_test == 0),
    n_masked_terminal_labels_for_ASR = asr$n_masked_terminal_labels,
    n_ASR_state0_total = feature_res$n_state0_total,
    n_ASR_state0_5_total = feature_res$n_state0_5_total,
    n_ASR_state1_total = feature_res$n_state1_total,
    n_ASR_state_other_total = feature_res$n_state_other_total,
    n_branches_screened_after_masking_and_heldout_removal = feature_res$n_branches_screened,
    n_screen_state0 = feature_res$n_screen_state0,
    n_screen_state1 = feature_res$n_screen_state1,
    n_heldout_terminal_branches_removed_from_ttest = feature_res$n_heldout_terminal_branches_removed_from_ttest,
    n_ttest_genes_tested = feature_res$n_ttest_genes_tested,
    n_fold_candidate_genes = length(fold_genes),
    stringsAsFactors = FALSE
  )

  null_return <- function(status, pred, dropped_all = NA_integer_, dropped_zero = NA_integer_, notes = "") {
    list(
      fold = cbind(
        base_fold,
        n_features_dropped_all_missing = dropped_all,
        n_features_dropped_zero_variance = dropped_zero,
        n_features_used_for_lasso = 0L,
        n_selected_predictors = 0L,
        lambda_min = NA_real_,
        fold_model_status = status,
        fold_seed = fold_seed,
        notes = notes,
        stringsAsFactors = FALSE
      ),
      oof = data.frame(
        trait = spec$trait,
        model = spec$model,
        frozen_model = spec$frozen_model,
        species = test_species$Species,
        genus = test_species$genus,
        trait_category = test_species$trait_category,
        marine_binary = test_species$marine_binary,
        binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA_real_, test_species$aquatic_v2),
        aquaticity_score = test_species$aquaticity_score,
        held_out_genus = fold_genus,
        foldwise_ASR_OOF_prediction = pred,
        foldwise_ASR_OOF_available = TRUE,
        exclusion_reason = "",
        fold_n_candidate_features_FDR = length(fold_genes),
        fold_n_features_after_preprocessing = 0L,
        fold_n_selected_predictors = 0L,
        fold_model_status = status,
        stringsAsFactors = FALSE
      ),
      selected = data.frame(),
      warnings = data.frame(
        trait = spec$trait, model = spec$model, fold_id = fold_id, held_out_genus = fold_genus,
        warning_type = "null_model",
        warning_message = status,
        stringsAsFactors = FALSE
      )
    )
  }

  if (length(unique(y_train)) < 2) {
    return(null_return(
      "null_intercept_only_training_one_class",
      rep(mean(y_train), length(y_test)),
      notes = "Training fold has only one response class."
    ))
  }
  if (length(fold_genes) == 0) {
    return(null_return(
      "null_intercept_only_no_FDR_features_in_training_foldwise_ASR_ttest",
      rep(mean(y_train), length(y_test)),
      dropped_all = 0L,
      dropped_zero = 0L,
      notes = "Fold-wise ASR + training-only t-test/FDR selected zero candidate genes."
    ))
  }

  x_train_raw <- t(gbi_mat_all[fold_genes, train_species$Branch, drop = FALSE])
  x_test_raw <- t(gbi_mat_all[fold_genes, test_species$Branch, drop = FALSE])
  train_means <- colMeans(x_train_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(train_means)
  dropped_all_missing <- sum(!keep_not_all_missing)
  x_train_raw <- x_train_raw[, keep_not_all_missing, drop = FALSE]
  x_test_raw <- x_test_raw[, keep_not_all_missing, drop = FALSE]
  train_means <- train_means[keep_not_all_missing]
  if (length(train_means) == 0) {
    return(null_return(
      "null_intercept_only_no_features_after_all_missing_drop",
      rep(mean(y_train), length(y_test)),
      dropped_all = dropped_all_missing,
      dropped_zero = 0L,
      notes = "No usable features remain after all-missing training drop."
    ))
  }
  for (k in seq_along(train_means)) {
    if (anyNA(x_train_raw[, k])) x_train_raw[is.na(x_train_raw[, k]), k] <- train_means[[k]]
    if (anyNA(x_test_raw[, k])) x_test_raw[is.na(x_test_raw[, k]), k] <- train_means[[k]]
  }
  scale_means <- colMeans(x_train_raw)
  scale_sds <- apply(x_train_raw, 2, stats::sd)
  keep_nonzero <- is.finite(scale_sds) & scale_sds > 0
  dropped_zero_variance <- sum(!keep_nonzero)
  x_train <- x_train_raw[, keep_nonzero, drop = FALSE]
  x_test <- x_test_raw[, keep_nonzero, drop = FALSE]
  scale_means <- scale_means[keep_nonzero]
  scale_sds <- scale_sds[keep_nonzero]
  if (ncol(x_train) == 0) {
    return(null_return(
      "null_intercept_only_no_features_after_zero_variance_drop",
      rep(mean(y_train), length(y_test)),
      dropped_all = dropped_all_missing,
      dropped_zero = dropped_zero_variance,
      notes = "No usable features remain after zero-variance training drop."
    ))
  }
  x_train <- sweep(sweep(x_train, 2, scale_means, "-"), 2, scale_sds, "/")
  x_test <- sweep(sweep(x_test, 2, scale_means, "-"), 2, scale_sds, "/")
  cv_res <- safe_cv_glmnet(as.matrix(x_train), y_train, fold_seed)
  lambda_min <- cv_res$cv$lambda.min
  fit <- glmnet(
    x = as.matrix(x_train),
    y = y_train,
    family = "binomial",
    lambda = lambda_min,
    standardize = FALSE
  )
  pred <- as.numeric(predict(fit, newx = as.matrix(x_test), type = "response")[, 1])
  beta <- as.matrix(fit$beta)[, 1]
  nonzero <- beta[beta != 0]
  selected <- if (length(nonzero) == 0) {
    data.frame()
  } else {
    data.frame(
      trait = spec$trait,
      model = spec$model,
      frozen_model = spec$frozen_model,
      fold_id = fold_id,
      held_out_genus = fold_genus,
      gene = names(nonzero),
      coefficient = as.numeric(nonzero),
      coefficient_sign = ifelse(nonzero > 0, "positive", "negative"),
      lambda_min = lambda_min,
      stringsAsFactors = FALSE
    )
  }
  warning_df <- if (length(cv_res$warnings) == 0) {
    data.frame()
  } else {
    data.frame(
      trait = spec$trait,
      model = spec$model,
      fold_id = fold_id,
      held_out_genus = fold_genus,
      warning_type = "cv.glmnet_warning",
      warning_message = cv_res$warnings,
      stringsAsFactors = FALSE
    )
  }
  status <- ifelse(length(nonzero) == 0, "fit_no_selected_predictors", "fit_nonzero_predictors")
  list(
    fold = cbind(
      base_fold,
      n_features_dropped_all_missing = dropped_all_missing,
      n_features_dropped_zero_variance = dropped_zero_variance,
      n_features_used_for_lasso = ncol(x_train),
      n_selected_predictors = length(nonzero),
      lambda_min = lambda_min,
      fold_model_status = status,
      fold_seed = fold_seed,
      notes = "Held-out genus terminal labels masked for ASR; held-out terminal branches removed from screening; fold-wise preprocessing; glmnet standardize=FALSE.",
      stringsAsFactors = FALSE
    ),
    oof = data.frame(
      trait = spec$trait,
      model = spec$model,
      frozen_model = spec$frozen_model,
      species = test_species$Species,
      genus = test_species$genus,
      trait_category = test_species$trait_category,
      marine_binary = test_species$marine_binary,
      binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA_real_, test_species$aquatic_v2),
      aquaticity_score = test_species$aquaticity_score,
      held_out_genus = fold_genus,
      foldwise_ASR_OOF_prediction = pred,
      foldwise_ASR_OOF_available = TRUE,
      exclusion_reason = "",
      fold_n_candidate_features_FDR = length(fold_genes),
      fold_n_features_after_preprocessing = ncol(x_train),
      fold_n_selected_predictors = length(nonzero),
      fold_model_status = status,
      stringsAsFactors = FALSE
    ),
    selected = selected,
    warnings = warning_df
  )
}

run_model <- function(spec, smoke = FALSE) {
  all_genera <- sort(unique(terminal$genus))
  if (smoke) {
    eligible_species <- terminal[spec$eligible, , drop = FALSE]
    choose_genus <- function(category = NULL, response = NULL, preferred = character()) {
      pool <- eligible_species
      if (!is.null(category)) pool <- pool[pool$trait_category == category, , drop = FALSE]
      if (!is.null(response)) pool <- pool[pool[[spec$response_col]] == response, , drop = FALSE]
      if (nrow(pool) == 0) return(NA_character_)
      for (p in preferred) if (p %in% pool$genus) return(p)
      sort(unique(pool$genus))[[1]]
    }
    chosen <- if (spec$trait == "marine_binary") {
      c(
        choose_genus(category = "marine", preferred = c("Orcinus", "Tursiops", "Phoca")),
        choose_genus(category = "non-marine aquatic", preferred = c("Hippopotamus", "Castor", "Ondatra", "Lutra")),
        choose_genus(category = "terrestrial", preferred = c("Acinonyx", "Ailuropoda")),
        choose_genus(category = "marine", preferred = c("Ursus", "Enhydra"))
      )
    } else {
      c(
        choose_genus(response = 1, category = "marine", preferred = c("Orcinus", "Tursiops")),
        choose_genus(response = 1, category = "non-marine aquatic", preferred = c("Hippopotamus", "Castor", "Ondatra", "Lutra")),
        choose_genus(response = 0, category = "terrestrial", preferred = c("Acinonyx", "Ailuropoda")),
        choose_genus(response = 1, category = "marine", preferred = c("Ursus", "Enhydra"))
      )
    }
    all_genera <- unique(chosen[!is.na(chosen)])
  }

  stamp("Running ", ifelse(smoke, "smoke ", ""), spec$model, " across ", length(all_genera), " genus folds")
  chunks <- split(all_genera, ceiling(seq_along(all_genera) / cores))
  results <- list()
  genus_index <- match(all_genera, sort(unique(terminal$genus)))
  for (chunk_i in seq_along(chunks)) {
    genera <- chunks[[chunk_i]]
    idx <- genus_index[match(genera, all_genera)]
    chunk_results <- parallel::mclapply(
      seq_along(genera),
      function(k) fit_fold(spec, genera[[k]], idx[[k]]),
      mc.cores = min(cores, length(genera))
    )
    results <- c(results, chunk_results)
    fold_df <- do.call(rbind, lapply(results, `[[`, "fold"))
    checkpoint <- ifelse(smoke, paste0(spec$trait, ".smoke.foldwise_ASR.fold_diagnostics.checkpoint.tsv"), paste0(spec$trait, ".foldwise_ASR.fold_diagnostics.checkpoint.tsv"))
    write_tsv(fold_df, file.path(work_dir, checkpoint))
    stamp(spec$model, ": completed ", length(results), "/", length(all_genera), " genus folds")
  }
  fold_diag <- do.call(rbind, lapply(results, `[[`, "fold"))
  oof <- do.call(rbind, lapply(results, `[[`, "oof"))
  selected <- do.call(rbind, lapply(results, `[[`, "selected"))
  warnings <- do.call(rbind, lapply(results, `[[`, "warnings"))
  if (is.null(selected)) selected <- data.frame()
  if (is.null(warnings)) warnings <- data.frame()

  if (!smoke) {
    excluded <- terminal[!spec$eligible, , drop = FALSE]
    if (nrow(excluded) > 0) {
      excluded_rows <- data.frame(
        trait = spec$trait,
        model = spec$model,
        frozen_model = spec$frozen_model,
        species = excluded$Species,
        genus = excluded$genus,
        trait_category = excluded$trait_category,
        marine_binary = excluded$marine_binary,
        binary_aquatic_endpoint = NA_real_,
        aquaticity_score = excluded$aquaticity_score,
        held_out_genus = excluded$genus,
        foldwise_ASR_OOF_prediction = NA_real_,
        foldwise_ASR_OOF_available = FALSE,
        exclusion_reason = "semi-aquatic aquatic_v2=0.5 excluded from binary aquatic-dependence LASSO/AUC",
        fold_n_candidate_features_FDR = NA_integer_,
        fold_n_features_after_preprocessing = NA_integer_,
        fold_n_selected_predictors = NA_integer_,
        fold_model_status = "excluded_no_binary_aquatic_endpoint",
        stringsAsFactors = FALSE
      )
      oof <- rbind(oof, excluded_rows)
    }
  }
  eval_oof <- oof[oof$foldwise_ASR_OOF_available, , drop = FALSE]
  y <- if (spec$trait == "marine_binary") eval_oof$marine_binary else eval_oof$binary_aquatic_endpoint
  auc <- auc_rank(y, eval_oof$foldwise_ASR_OOF_prediction)
  list(spec = spec, fold_diag = fold_diag, oof = oof, selected = selected, warnings = warnings, auc = auc)
}

if (mode == "smoke") {
  stamp("Starting fold-wise ASR smoke test")
  smoke_results <- lapply(model_specs, run_model, smoke = TRUE)
  smoke_folds <- do.call(rbind, lapply(smoke_results, `[[`, "fold_diag"))
  smoke_oof <- do.call(rbind, lapply(smoke_results, `[[`, "oof"))
  smoke_warnings <- do.call(rbind, lapply(smoke_results, `[[`, "warnings"))
  write_tsv(smoke_folds, file.path(qc_dir, "foldwise_ASR_smoke_fold_level_summary.tsv"))
  write_tsv(smoke_oof, file.path(qc_dir, "foldwise_ASR_smoke_OOF_predictions.tsv"))
  if (is.null(smoke_warnings) || nrow(smoke_warnings) == 0) {
    smoke_warnings <- data.frame(trait = character(), model = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())
  }
  write_tsv(smoke_warnings, file.path(qc_dir, "foldwise_ASR_smoke_warning_log.tsv"))
  stamp("Fold-wise ASR smoke test complete")
  quit(save = "no")
}

stamp("Starting full fold-wise ASR sensitivity run with cores=", cores)
results <- lapply(model_specs, run_model, smoke = FALSE)
fold_diagnostics <- do.call(rbind, lapply(results, `[[`, "fold_diag"))
all_oof <- do.call(rbind, lapply(results, `[[`, "oof"))
selected_by_fold <- do.call(rbind, lapply(results, `[[`, "selected"))
warning_log <- do.call(rbind, lapply(results, `[[`, "warnings"))
if (is.null(selected_by_fold)) selected_by_fold <- data.frame()
if (is.null(warning_log) || nrow(warning_log) == 0) {
  warning_log <- data.frame(trait = character(), model = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())
}

write_tsv(fold_diagnostics, file.path(work_dir, "foldwise_ASR_fold_level_summary.tsv"))
write_tsv(all_oof, file.path(work_dir, "foldwise_ASR_OOF_predictions.tsv"))
write_tsv(selected_by_fold, file.path(work_dir, "foldwise_ASR_selected_predictors_by_fold.tsv"))
write_tsv(warning_log, file.path(work_dir, "foldwise_ASR_warning_log.tsv"))

frozen_selected_counts <- if (nrow(phase12b_selected) > 0) {
  aggregate(gene ~ trait + model + fold_id + held_out_genus, data = phase12b_selected, FUN = length)
} else {
  data.frame(trait = character(), model = character(), fold_id = integer(), held_out_genus = character(), gene = integer())
}
names(frozen_selected_counts)[names(frozen_selected_counts) == "gene"] <- "frozen_n_selected_predictors"

fold_compare <- merge(
  fold_diagnostics,
  phase12b_fold[, c("trait", "model", "fold_id", "held_out_genus", "n_fold_candidate_genes", "n_selected_predictors")],
  by.x = c("trait", "frozen_model", "fold_id", "held_out_genus"),
  by.y = c("trait", "model", "fold_id", "held_out_genus"),
  all.x = TRUE,
  suffixes = c("_foldwise_ASR", "_frozen_ASR")
)
names(fold_compare)[names(fold_compare) == "n_fold_candidate_genes_frozen_ASR"] <- "frozen_ASR_n_candidate_genes"
names(fold_compare)[names(fold_compare) == "n_selected_predictors_frozen_ASR"] <- "frozen_ASR_n_selected_predictors"
fold_compare$delta_candidate_genes_foldwise_minus_frozen <- fold_compare$n_fold_candidate_genes - fold_compare$frozen_ASR_n_candidate_genes
fold_compare$delta_selected_predictors_foldwise_minus_frozen <- fold_compare$n_selected_predictors - fold_compare$frozen_ASR_n_selected_predictors
write_tsv(fold_compare, file.path(work_dir, "foldwise_ASR_fold_level_summary.tsv"))

auc_summary <- do.call(rbind, lapply(results, function(res) {
  spec <- res$spec
  fold <- res$fold_diag
  eval_oof <- res$oof[res$oof$foldwise_ASR_OOF_available, , drop = FALSE]
  y <- if (spec$trait == "marine_binary") eval_oof$marine_binary else eval_oof$binary_aquatic_endpoint
  frozen_auc_row <- phase12b_auc[phase12b_auc$model == spec$frozen_model, , drop = FALSE]
  data.frame(
    trait = spec$trait,
    trait_axis = spec$trait_axis,
    foldwise_ASR_model = spec$model,
    frozen_ASR_model = spec$frozen_model,
    n_species_evaluated = sum(res$oof$foldwise_ASR_OOF_available),
    n_positive = sum(y == 1, na.rm = TRUE),
    n_negative = sum(y == 0, na.rm = TRUE),
    semi_aquatic_handling = spec$semi_aquatic_handling,
    AUC_frozen_ASR_nested_ttest_raw = frozen_auc_row$AUC_nested_ttest[[1]],
    AUC_frozen_ASR_nested_ttest_public = spec$frozen_auc_public,
    AUC_foldwise_ASR_nested_ttest = res$auc,
    delta_AUC_foldwise_minus_frozen = res$auc - frozen_auc_row$AUC_nested_ttest[[1]],
    n_LOOCV_folds = sum(fold$n_test_samples > 0),
    n_folds_zero_candidate_features = sum(fold$n_test_samples > 0 & fold$n_fold_candidate_genes == 0, na.rm = TRUE),
    n_folds_null_intercept_only = sum(fold$n_test_samples > 0 & grepl("^null_intercept_only", fold$fold_model_status)),
    median_candidate_features_per_fold = stats::median(fold$n_fold_candidate_genes[fold$n_test_samples > 0], na.rm = TRUE),
    IQR_candidate_features_per_fold = iqr_text(fold$n_fold_candidate_genes[fold$n_test_samples > 0], digits = 0),
    median_selected_predictors_per_fold = stats::median(fold$n_selected_predictors[fold$n_test_samples > 0], na.rm = TRUE),
    IQR_selected_predictors_per_fold = iqr_text(fold$n_selected_predictors[fold$n_test_samples > 0], digits = 0),
    evidence_level = "fold-wise ASR sensitivity: held-out genus terminal trait labels masked before deterministic maximum-parsimony branch-state reconstruction; supervised t-test/FDR and terminal-only LASSO nested inside genus folds",
    material_change_flag = ifelse(abs(res$auc - frozen_auc_row$AUC_nested_ttest[[1]]) > 0.02, "REVIEW_MATERIAL_AUC_CHANGE", "PASS_NO_MATERIAL_AUC_CHANGE"),
    stringsAsFactors = FALSE
  )
}))
write_tsv(auc_summary, file.path(work_dir, "foldwise_ASR_vs_frozen_ASR_AUC_summary.tsv"))

candidate_summary <- do.call(rbind, lapply(split(fold_compare[fold_compare$n_test_samples > 0, , drop = FALSE], fold_compare$trait[fold_compare$n_test_samples > 0]), function(df) {
  data.frame(
    trait = unique(df$trait),
    foldwise_ASR_model = unique(df$model),
    frozen_ASR_model = unique(df$frozen_model),
    n_evaluable_folds = nrow(df),
    frozen_median_candidate_genes = stats::median(df$frozen_ASR_n_candidate_genes, na.rm = TRUE),
    frozen_IQR_candidate_genes = iqr_text(df$frozen_ASR_n_candidate_genes, digits = 0),
    foldwise_median_candidate_genes = stats::median(df$n_fold_candidate_genes, na.rm = TRUE),
    foldwise_IQR_candidate_genes = iqr_text(df$n_fold_candidate_genes, digits = 0),
    median_delta_candidate_genes = stats::median(df$delta_candidate_genes_foldwise_minus_frozen, na.rm = TRUE),
    n_folds_candidate_count_changed = sum(df$delta_candidate_genes_foldwise_minus_frozen != 0, na.rm = TRUE),
    max_abs_delta_candidate_genes = max(abs(df$delta_candidate_genes_foldwise_minus_frozen), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv(candidate_summary, file.path(work_dir, "foldwise_ASR_candidate_gene_count_summary.tsv"))

selected_summary <- do.call(rbind, lapply(split(fold_compare[fold_compare$n_test_samples > 0, , drop = FALSE], fold_compare$trait[fold_compare$n_test_samples > 0]), function(df) {
  data.frame(
    trait = unique(df$trait),
    foldwise_ASR_model = unique(df$model),
    frozen_ASR_model = unique(df$frozen_model),
    n_evaluable_folds = nrow(df),
    frozen_median_selected_predictors = stats::median(df$frozen_ASR_n_selected_predictors, na.rm = TRUE),
    frozen_IQR_selected_predictors = iqr_text(df$frozen_ASR_n_selected_predictors, digits = 0),
    foldwise_median_selected_predictors = stats::median(df$n_selected_predictors, na.rm = TRUE),
    foldwise_IQR_selected_predictors = iqr_text(df$n_selected_predictors, digits = 0),
    median_delta_selected_predictors = stats::median(df$delta_selected_predictors_foldwise_minus_frozen, na.rm = TRUE),
    n_folds_selected_count_changed = sum(df$delta_selected_predictors_foldwise_minus_frozen != 0, na.rm = TRUE),
    max_abs_delta_selected_predictors = max(abs(df$delta_selected_predictors_foldwise_minus_frozen), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write_tsv(selected_summary, file.path(work_dir, "foldwise_ASR_selected_predictor_count_summary.tsv"))

methods_sentence <- if (all(auc_summary$material_change_flag == "PASS_NO_MATERIAL_AUC_CHANGE")) {
  sprintf(
    "As a sensitivity check, deterministic branch-state reconstruction was repeated inside each outer genus-level fold after masking held-out terminal labels; baseline AUCs were similar to the frozen-ASR nested analysis (marine %.3f vs %.3f; binary aquatic-dependence %.3f vs %.3f), indicating that fixed ASR labels did not materially affect the baseline gLOOCV results.",
    auc_summary$AUC_foldwise_ASR_nested_ttest[auc_summary$trait == "marine_binary"],
    auc_summary$AUC_frozen_ASR_nested_ttest_raw[auc_summary$trait == "marine_binary"],
    auc_summary$AUC_foldwise_ASR_nested_ttest[auc_summary$trait == "binary_aquatic_dependence"],
    auc_summary$AUC_frozen_ASR_nested_ttest_raw[auc_summary$trait == "binary_aquatic_dependence"]
  )
} else {
  "Fold-wise ASR changed at least one baseline AUC by more than 0.02 relative to the frozen-ASR nested analysis; main baseline AUCs should be reviewed before manuscript integration."
}

notes <- c(
  "# Fold-Wise ASR Missing-Label Handling",
  "",
  "This sensitivity check uses `castor::hsp_max_parsimony`, the castor maximum-parsimony routine for hidden/unknown tip states.",
  "",
  "For each outer genus-level fold:",
  "",
  "1. Terminal trait labels for the held-out genus were set to `NA` for ASR only.",
  "2. Deterministic maximum-parsimony reconstruction used the same support tree, all-equal transition costs, edge exponent 0 and first-state tie handling as the frozen-ASR pipeline.",
  "3. Known terminal labels were preserved directly in the branch table, matching the Stage 1 branch-state mapping convention.",
  "4. Node states were mapped to branches using the reconstructed/predicted state of the descendant node of each branch.",
  "5. Held-out genus terminal branches were then removed from the branch-level t-test/FDR screening table.",
  "6. Branches with intermediate `0.5` states were excluded from binary branch-level contrasts, matching the Stage 1 convention.",
  "",
  "A no-missing validation confirmed that `hsp_max_parsimony` reproduces the existing frozen deterministic ASR branch-state vectors exactly for the marine and binary aquatic-dependence baseline traits.",
  "",
  "Suggested note if accepted:",
  "",
  methods_sentence,
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(notes, file.path(work_dir, "notes_on_ASR_missing_label_handling.md"))

readme <- c(
  "# Fold-Wise ASR Sensitivity Check",
  "",
  "Purpose: evaluate whether keeping deterministic ASR branch-state labels frozen across folds materially affects the main nested gLOOCV baseline AUCs.",
  "",
  "No GBI, SplitAligner, enrichment or manuscript files were changed. This run reuses the endpoint-fix GBI matrix and reruns only deterministic ASR branch-state assignment inside each outer fold as a sensitivity check.",
  "",
  "Required outputs:",
  "",
  "- `foldwise_ASR_vs_frozen_ASR_AUC_summary.tsv`",
  "- `foldwise_ASR_fold_level_summary.tsv`",
  "- `foldwise_ASR_candidate_gene_count_summary.tsv`",
  "- `foldwise_ASR_selected_predictor_count_summary.tsv`",
  "- `notes_on_ASR_missing_label_handling.md`",
  "",
  "Additional audit outputs:",
  "",
  "- `foldwise_ASR_OOF_predictions.tsv`",
  "- `foldwise_ASR_selected_predictors_by_fold.tsv`",
  "- `foldwise_ASR_warning_log.tsv`",
  "- `qc/foldwise_ASR_hsp_vs_frozen_validation.tsv`",
  "",
  "Bottom line:",
  "",
  methods_sentence,
  "",
  paste0("Runtime: ", sprintf("%.1f minutes", as.numeric(difftime(Sys.time(), script_start, units = "mins")))),
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(readme, file.path(work_dir, "README_foldwise_ASR_sensitivity.md"))

qc <- data.frame(
  check = c(
    "hsp_no_missing_reproduces_frozen_ASR_marine",
    "hsp_no_missing_reproduces_frozen_ASR_binary_aquatic",
    "held_out_terminal_labels_masked_for_ASR",
    "held_out_terminal_branches_removed_from_ttest",
    "foldwise_supervised_ttest_FDR",
    "foldwise_training_only_imputation_scaling_lambda_fit",
    "glmnet_standardize_FALSE",
    "AUC_rank_based_pooled_OOF",
    "material_AUC_change_threshold_0.02"
  ),
  status = c(
    ifelse(hsp_validation$status[hsp_validation$trait == "marine_binary"] == "PASS", "PASS", "FAIL"),
    ifelse(hsp_validation$status[hsp_validation$trait == "binary_aquatic_dependence"] == "PASS", "PASS", "FAIL"),
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    ifelse(all(auc_summary$material_change_flag == "PASS_NO_MATERIAL_AUC_CHANGE"), "PASS", "REVIEW")
  ),
  notes = c(
    "No-missing hsp_max_parsimony branch-state vector matched frozen Stage 1 ASR exactly for marine.",
    "No-missing hsp_max_parsimony branch-state vector matched frozen Stage 1 ASR exactly for binary aquatic-dependence.",
    "For each fold, held-out genus terminal labels are NA in hsp_max_parsimony only.",
    "After fold-wise ASR, held-out genus terminal branches are removed from branch-level screening.",
    "Welch t-test/FDR candidate selection is repeated inside each fold.",
    "LASSO preprocessing and cv.glmnet/glmnet use training terminal samples only.",
    "Manual scaling is used; glmnet/cv.glmnet standardize=FALSE.",
    "AUC is the same rank-based pooled out-of-fold statistic.",
    "AUC changes <= 0.02 are treated as no material baseline change for this sensitivity check."
  ),
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(qc_dir, "foldwise_ASR_sensitivity_QC.tsv"))

stamp("Fold-wise ASR sensitivity complete")
