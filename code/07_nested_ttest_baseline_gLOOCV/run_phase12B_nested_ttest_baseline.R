#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(gridExtra)
  library(parallel)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness")
code_dir <- file.path(work_dir, "code")
design_dir <- file.path(work_dir, "design_audit")
smoke_dir <- file.path(work_dir, "smoke_test")
smoke_log_dir <- file.path(smoke_dir, "nested_ttest_smoke_test_logs")
out_dir <- file.path(work_dir, "nested_ttest_baseline")
qc_dir <- file.path(work_dir, "qc")
decision_dir <- file.path(work_dir, "decision")
fig_dir <- file.path(work_dir, "figures")
dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(smoke_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(smoke_log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(decision_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}
mode <- get_arg("--mode", "full")
cores <- as.integer(get_arg("--cores", min(10L, max(1L, parallel::detectCores() - 2L))))
if (!is.finite(cores) || cores < 1L) cores <- 1L
seed_base <- as.integer(get_arg("--seed", 20260523L))
alpha <- as.numeric(get_arg("--alpha", 0.01))
if (!mode %in% c("smoke", "full")) stop("--mode must be smoke or full")
script_start <- Sys.time()

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  asr_branch_states = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "ancestral_state_plots", "endpointfix_stage1_ancestral_state_branch_table.tsv"),
  global_marine_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_marine_binary", "marine_binary.mammal.FDR0.01.n1559.t_test.txt"),
  global_aquatic_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_aquatic_v2", "aquatic_v2.mammal.FDR0.01.n1227.t_test.txt"),
  phase11_auc = file.path(root_dir, "10_reviewer_risk_controls", "11_corrected_preprocessing_terminal_LASSO", "baseline_corrected", "corrected_preprocessing_AUC_summary.tsv"),
  phase11_oof = file.path(root_dir, "10_reviewer_risk_controls", "11_corrected_preprocessing_terminal_LASSO", "baseline_corrected", "corrected_preprocessing_OOF_predictions_all_species.tsv"),
  phase11_category = file.path(root_dir, "10_reviewer_risk_controls", "11_corrected_preprocessing_terminal_LASSO", "baseline_corrected", "legacy_vs_corrected_category_OOF_summary.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

write_tsv <- function(x, file) {
  utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}
read_tsv <- function(file, ...) {
  utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
}
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

roc_points <- function(y, p, model, auc_value) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]
  p <- p[keep]
  thresholds <- c(Inf, sort(unique(p), decreasing = TRUE), -Inf)
  out <- lapply(thresholds, function(th) {
    pred_pos <- p >= th
    tp <- sum(pred_pos & y == 1)
    fp <- sum(pred_pos & y == 0)
    fn <- sum(!pred_pos & y == 1)
    tn <- sum(!pred_pos & y == 0)
    data.frame(
      model = model,
      threshold = th,
      FPR = ifelse((fp + tn) > 0, fp / (fp + tn), NA_real_),
      TPR = ifelse((tp + fn) > 0, tp / (tp + fn), NA_real_),
      AUC = auc_value,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, out)
}

iqr_text <- function(x, digits = 0) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  fmt <- paste0("%.", digits, "f-%.", digits, "f")
  sprintf(fmt, qs[[1]], qs[[2]])
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

subtract_stats <- function(base, heldout) {
  list(
    n = base$n - heldout$n,
    sum = base$sum - heldout$sum,
    sumsq = base$sumsq - heldout$sumsq
  )
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
asr_states <- read_tsv(paths$asr_branch_states)
phase11_auc <- read_tsv(paths$phase11_auc)
phase11_oof <- read_tsv(paths$phase11_oof)

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

global_genes <- list(
  marine_binary = read_tsv(paths$global_marine_features)$gene,
  binary_aquatic_dependence = read_tsv(paths$global_aquatic_features)$gene
)

model_specs <- list(
  list(
    trait = "marine_binary",
    trait_axis = "marine specialization",
    model = "marine_binary_nested_ttest_baseline",
    phase11_model = "marine_binary_corrected_preprocessing_gLOOCV",
    run_id = "fix_marine_binary",
    trait_column = "marine_binary",
    response_col = "marine_binary",
    eligible = rep(TRUE, nrow(terminal)),
    global_candidate_genes = global_genes$marine_binary,
    semi_aquatic_handling = "retained as non-marine/background 0 for marine binary model"
  ),
  list(
    trait = "binary_aquatic_dependence",
    trait_axis = "binary aquatic dependence",
    model = "binary_aquatic_dependence_nested_ttest_baseline",
    phase11_model = "binary_aquatic_dependence_corrected_preprocessing_gLOOCV",
    run_id = "fix_aquatic_v2",
    trait_column = "aquatic_v2",
    response_col = "aquatic_v2",
    eligible = terminal$aquatic_v2 != 0.5,
    global_candidate_genes = global_genes$binary_aquatic_dependence,
    semi_aquatic_handling = "semi-aquatic aquatic_v2=0.5 labels excluded from training, testing, imputation, scaling, and AUC"
  )
)

make_trait_cache <- function(spec) {
  states <- asr_states[asr_states$run_id == spec$run_id & asr_states$trait_column == spec$trait_column, , drop = FALSE]
  if (nrow(states) == 0) stop("No ASR branch-state rows for ", spec$run_id)
  states <- states[match(branch$Branch, states$old_branch_label), , drop = FALSE]
  if (!all(states$old_branch_label == branch$Branch)) stop("ASR branch-state order mismatch for ", spec$run_id)
  state <- states$state
  cols0 <- branch$Branch[state == 0]
  cols1 <- branch$Branch[state == 1]
  cols05 <- branch$Branch[state == 0.5]
  base0 <- group_stats(gbi_mat_all, cols0)
  base1 <- group_stats(gbi_mat_all, cols1)
  terminal_state <- setNames(state[match(terminal$Branch, branch$Branch)], terminal$Species)
  global_t <- welch_from_stats(base0, base1)
  global_tested <- global_t[is.finite(global_t$pvalue), , drop = FALSE]
  global_idx <- bh_old_while(global_tested$pvalue, alpha = alpha)
  global_selected <- global_tested$gene[global_idx]
  list(
    states = states,
    state = state,
    cols0 = cols0,
    cols1 = cols1,
    cols05 = cols05,
    base0 = base0,
    base1 = base1,
    terminal_state = terminal_state,
    global_validation = data.frame(
      trait = spec$trait,
      n_global_ttest_genes_tested_vectorized = nrow(global_tested),
      n_global_FDR_candidates_vectorized = length(global_selected),
      n_global_FDR_candidates_archived = length(spec$global_candidate_genes),
      candidate_set_matches_archived = setequal(global_selected, spec$global_candidate_genes),
      notes = "Vectorized Welch t-test using frozen deterministic ASR branch states and old BH while-loop rule.",
      stringsAsFactors = FALSE
    )
  )
}

trait_caches <- lapply(model_specs, make_trait_cache)
names(trait_caches) <- vapply(model_specs, `[[`, character(1), "trait")

design_rows <- do.call(rbind, lapply(model_specs, function(spec) {
  cache <- trait_caches[[spec$trait]]
  data.frame(
    model = spec$model,
    trait_axis = ifelse(spec$trait == "marine_binary", "marine specialization", "binary aquatic dependence"),
    frozen_GBI_matrix = basename(paths$gbi),
    frozen_branch_state_source = basename(paths$asr_branch_states),
    training_data_only_definition_for_branch_ttest = "all frozen branch-state-labelled branches except terminal branches belonging to the held-out genus",
    held_out_genus_terminal_branches_excluded_from_ttest = "yes",
    internal_branches_included_in_ttest = "yes",
    internal_branch_justification = "retained to match the original deterministic branch-level discovery design; internal branches are frozen ASR-labelled branch observations and are never used as LASSO samples",
    branch_state_ASR_labels_frozen_globally = "yes",
    GBI_values_frozen_globally = "yes",
    same_t_statistic_direction_as_Stage1 = "yes; tvalue = mean(state0) - mean(state1)",
    same_BH_FDR_rule_as_Stage1 = "yes; old while-loop rule p[k] < k/m*alpha after sorting p-values",
    semi_aquatic_handling = spec$semi_aquatic_handling,
    n_total_branches = length(cache$state),
    n_terminal_branches = nrow(terminal),
    n_internal_branches = sum(branch$Species == "internal"),
    n_state0_branches_global = sum(cache$state == 0),
    n_state1_branches_global = sum(cache$state == 1),
    n_state0_5_branches_global = sum(cache$state == 0.5),
    status = "PASS_DESIGN_DEFINED",
    stringsAsFactors = FALSE
  )
}))
write_tsv(design_rows, file.path(design_dir, "nested_ttest_training_branch_definition.tsv"))

design_md <- c(
  "# Phase 12B Nested Supervised T-Test Design Audit",
  "",
  "This audit is written before full baseline execution.",
  "",
  "## Frozen Inputs",
  "",
  "- Endpoint-fix GBI matrix is frozen.",
  "- Old-label branch coordinates are frozen.",
  "- Deterministic branch-state/ASR table from the current global t-test pipeline is frozen.",
  "- Terminal-only LASSO sample definition and Phase 11 corrected preprocessing are retained.",
  "",
  "## Changed Step",
  "",
  "Only supervised branch-level t-test/FDR candidate-gene selection is nested inside each outer genus-level LOOCV fold.",
  "",
  "## Answers to Control-Chat Questions",
  "",
  "**How is training data only defined for branch-level t-test?**",
  "",
  "For each outer held-out genus fold, training branch data are the frozen ASR-labelled branch observations after excluding terminal branches belonging to the held-out genus. Internal branches are retained as frozen branch-level discovery observations.",
  "",
  "**Are held-out genus terminal branches excluded from the t-test?**",
  "",
  "Yes. Their terminal branch columns are subtracted from the fold-level t-test/FDR discovery statistics before candidate genes are selected.",
  "",
  "**Are internal branches included or excluded?**",
  "",
  "Internal branches are included in the fold-level t-test/FDR discovery screen, matching the original branch-level discovery design. They are not used in imputation, scaling, glmnet fitting, lambda selection, prediction, or AUC.",
  "",
  "**If internal branches are included, what is the justification?**",
  "",
  "The accepted project design treats the branch-level t-test as deterministic ASR-based discovery over frozen branch states. The reviewer-risk issue is supervised candidate selection seeing held-out terminal genus information; therefore the robustness run removes held-out terminal branches while retaining the original frozen internal branch-level discovery layer.",
  "",
  "**Are branch-state/ASR labels frozen globally?**",
  "",
  "Yes. No ASR rerun is performed in Phase 12B.",
  "",
  "**Are GBI values frozen globally?**",
  "",
  "Yes. No GBI, SplitAligner, or branch matrix step is rerun.",
  "",
  "**Does the nested t-test use the same t-statistic direction and BH/FDR rule as Stage 1?**",
  "",
  "Yes. The t statistic is `mean(state0) - mean(state1)`, matching Stage 1. FDR candidate selection uses the old Stage 1 sorted-p-value while-loop rule `p[k] < k/m*alpha` at alpha = 0.01.",
  "",
  "**How are semi-aquatic labels handled for the binary aquatic-dependence model?**",
  "",
  "Semi-aquatic aquatic_v2=0.5 labels are excluded from binary aquatic LASSO training, testing, imputation, scaling, and AUC. Branches with state 0.5 are also excluded from the branch-level t-test screen, matching Stage 1.",
  "",
  "## Execution Gate",
  "",
  "The internal-branch rule is explicit and follows the current accepted design, so full execution may proceed if the smoke test passes.",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(design_md, file.path(design_dir, "nested_ttest_design_audit.md"))

select_fold_features <- function(spec, cache, held_out_genus) {
  held <- terminal[terminal$genus == held_out_genus, , drop = FALSE]
  held_cols0 <- held$Branch[cache$terminal_state[held$Species] == 0]
  held_cols1 <- held$Branch[cache$terminal_state[held$Species] == 1]
  h0 <- group_stats(gbi_mat_all, held_cols0)
  h1 <- group_stats(gbi_mat_all, held_cols1)
  g0 <- subtract_stats(cache$base0, h0)
  g1 <- subtract_stats(cache$base1, h1)
  tt <- welch_from_stats(g0, g1)
  tt_tested <- tt[is.finite(tt$pvalue), , drop = FALSE]
  tt_tested <- tt_tested[order(tt_tested$pvalue, decreasing = FALSE), , drop = FALSE]
  idx <- bh_old_while(tt_tested$pvalue, alpha = alpha)
  genes <- tt_tested$gene[idx]
  list(
    genes = genes,
    ttest = tt_tested,
    n_branches_screened = sum((cache$state == 0 | cache$state == 1) & !(branch$Branch %in% c(held_cols0, held_cols1))),
    n_heldout_terminal_branches_removed_from_ttest = length(held_cols0) + length(held_cols1),
    n_ttest_genes_tested = nrow(tt_tested)
  )
}

fit_fold <- function(spec, cache, fold_genus, fold_index) {
  fold_seed <- seed_base + fold_index + ifelse(spec$trait == "marine_binary", 100000L, 200000L)
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  test_idx <- which(eligible_species$genus == fold_genus)
  if (length(test_idx) == 0) {
    return(list(
      fold = data.frame(
        trait = spec$trait, model = spec$model, fold_id = fold_index, held_out_genus = fold_genus,
        n_training_samples = nrow(eligible_species), n_test_samples = 0,
        n_train_positive = sum(eligible_species[[spec$response_col]] == 1),
        n_train_negative = sum(eligible_species[[spec$response_col]] == 0),
        n_test_positive = 0, n_test_negative = 0,
        n_ttest_genes_tested = NA_integer_,
        n_fold_candidate_genes = 0L,
        n_overlap_with_global_FDR_candidates = 0L,
        prop_nested_candidates_in_global_FDR = NA_real_,
        prop_global_FDR_recovered_in_nested_candidates = 0,
        n_heldout_terminal_branches_removed_from_ttest = 0L,
        n_branches_screened_in_training_ttest = NA_integer_,
        n_features_dropped_all_missing = NA_integer_,
        n_features_dropped_zero_variance = NA_integer_,
        n_features_used_for_lasso = 0L,
        n_selected_predictors = 0L,
        lambda_min = NA_real_,
        fold_model_status = "no_evaluable_test_species_after_exclusion",
        corrected_preprocessing_boundary = "not_run_no_evaluable_test_species",
        fold_seed = fold_seed,
        notes = "No evaluable held-out species for this model after exclusion.",
        stringsAsFactors = FALSE
      ),
      oof = data.frame(),
      selected = data.frame(),
      warnings = data.frame(
        trait = spec$trait, model = spec$model, fold_id = fold_index, held_out_genus = fold_genus,
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
  if (!all(y_train %in% c(0, 1)) || !all(y_test %in% c(0, 1))) {
    stop("Non-binary response in ", spec$model, " fold ", fold_genus)
  }

  feature_res <- select_fold_features(spec, cache, fold_genus)
  fold_genes <- feature_res$genes
  overlap <- length(intersect(fold_genes, spec$global_candidate_genes))
  base_fold <- data.frame(
    trait = spec$trait, model = spec$model, fold_id = fold_index, held_out_genus = fold_genus,
    n_training_samples = length(y_train), n_test_samples = length(y_test),
    n_train_positive = sum(y_train == 1), n_train_negative = sum(y_train == 0),
    n_test_positive = sum(y_test == 1), n_test_negative = sum(y_test == 0),
    n_ttest_genes_tested = feature_res$n_ttest_genes_tested,
    n_fold_candidate_genes = length(fold_genes),
    n_overlap_with_global_FDR_candidates = overlap,
    prop_nested_candidates_in_global_FDR = ifelse(length(fold_genes) > 0, overlap / length(fold_genes), NA_real_),
    prop_global_FDR_recovered_in_nested_candidates = overlap / length(spec$global_candidate_genes),
    n_heldout_terminal_branches_removed_from_ttest = feature_res$n_heldout_terminal_branches_removed_from_ttest,
    n_branches_screened_in_training_ttest = feature_res$n_branches_screened,
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
        corrected_preprocessing_boundary = "training-only imputation/scaling/lambda/model fitting; terminal species only",
        fold_seed = fold_seed,
        notes = notes,
        stringsAsFactors = FALSE
      ),
      oof = data.frame(
        trait = spec$trait, model = spec$model,
        species = test_species$Species, genus = test_species$genus,
        trait_category = test_species$trait_category,
        marine_binary = test_species$marine_binary,
        binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2),
        aquaticity_score = test_species$aquaticity_score,
        held_out_genus = fold_genus,
        globalFDR_corrected_OOF_prediction = NA_real_,
        nested_ttest_OOF_prediction = pred,
        prediction_delta_nested_minus_globalFDR = NA_real_,
        nested_OOF_available = TRUE,
        exclusion_reason = "",
        fold_n_ttest_genes_tested = feature_res$n_ttest_genes_tested,
        fold_n_candidate_genes_from_training_ttest = length(fold_genes),
        fold_n_features_dropped_all_missing = dropped_all,
        fold_n_features_dropped_zero_variance = dropped_zero,
        fold_n_features_used = 0L,
        fold_n_selected_predictors = 0L,
        fold_lambda = NA_real_,
        fold_model_status = status,
        stringsAsFactors = FALSE
      ),
      selected = data.frame(),
      warnings = data.frame(
        trait = spec$trait, model = spec$model, fold_id = fold_index, held_out_genus = fold_genus,
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
      "null_intercept_only_no_FDR_features_in_training_ttest",
      rep(mean(y_train), length(y_test)),
      dropped_all = 0L,
      dropped_zero = 0L,
      notes = "Training-only t-test/FDR selected zero candidate genes."
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
      fold_id = fold_index,
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
      trait = spec$trait, model = spec$model, fold_id = fold_index, held_out_genus = fold_genus,
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
      corrected_preprocessing_boundary = "training-only imputation/scaling/lambda/model fitting; terminal species only",
      fold_seed = fold_seed,
      notes = "Fold-specific training t-test/FDR features; fold-wise preprocessing; glmnet standardize=FALSE.",
      stringsAsFactors = FALSE
    ),
    oof = data.frame(
      trait = spec$trait, model = spec$model,
      species = test_species$Species, genus = test_species$genus,
      trait_category = test_species$trait_category,
      marine_binary = test_species$marine_binary,
      binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2),
      aquaticity_score = test_species$aquaticity_score,
      held_out_genus = fold_genus,
      globalFDR_corrected_OOF_prediction = NA_real_,
      nested_ttest_OOF_prediction = pred,
      prediction_delta_nested_minus_globalFDR = NA_real_,
      nested_OOF_available = TRUE,
      exclusion_reason = "",
      fold_n_ttest_genes_tested = feature_res$n_ttest_genes_tested,
      fold_n_candidate_genes_from_training_ttest = length(fold_genes),
      fold_n_features_dropped_all_missing = dropped_all_missing,
      fold_n_features_dropped_zero_variance = dropped_zero_variance,
      fold_n_features_used = ncol(x_train),
      fold_n_selected_predictors = length(nonzero),
      fold_lambda = lambda_min,
      fold_model_status = status,
      stringsAsFactors = FALSE
    ),
    selected = selected,
    warnings = warning_df
  )
}

run_model <- function(spec, cache, smoke = FALSE) {
  all_genera <- sort(unique(terminal$genus))
  if (smoke) {
    eligible_species <- terminal[spec$eligible, , drop = FALSE]
    choose_genus <- function(category = NULL, response = NULL, preferred = character()) {
      pool <- eligible_species
      if (!is.null(category)) pool <- pool[pool$trait_category == category, , drop = FALSE]
      if (!is.null(response)) pool <- pool[pool[[spec$response_col]] == response, , drop = FALSE]
      if (nrow(pool) == 0) return(NA_character_)
      for (p in preferred) {
        if (p %in% pool$genus) return(p)
      }
      sort(unique(pool$genus))[[1]]
    }
    chosen <- if (spec$trait == "marine_binary") {
      c(
        choose_genus(category = "marine", preferred = c("Orcinus", "Tursiops", "Phoca")),
        choose_genus(category = "non-marine aquatic", preferred = c("Hippopotamus", "Castor", "Ondatra", "Lutra")),
        choose_genus(category = "terrestrial", preferred = c("Acinonyx", "Ailuropoda")),
        choose_genus(category = "marine", preferred = c("Ursus", "Enhydra")),
        choose_genus(category = "semi-aquatic", preferred = c("Tapirus", "Alces"))
      )
    } else {
      c(
        choose_genus(response = 1, category = "marine", preferred = c("Orcinus", "Tursiops")),
        choose_genus(response = 1, category = "non-marine aquatic", preferred = c("Hippopotamus", "Castor", "Ondatra", "Lutra")),
        choose_genus(response = 0, category = "terrestrial", preferred = c("Acinonyx", "Ailuropoda")),
        choose_genus(response = 1, category = "marine", preferred = c("Ursus", "Enhydra")),
        choose_genus(response = 0, category = "terrestrial", preferred = c("Hydrochoerus", "Cavia"))
      )
    }
    chosen <- unique(chosen[!is.na(chosen)])
    all_genera <- all_genera[all_genera %in% chosen]
  }
  stamp("Running ", ifelse(smoke, "smoke ", ""), "nested-t-test baseline for ", spec$model, " across ", length(all_genera), " genus folds")
  chunks <- split(all_genera, ceiling(seq_along(all_genera) / cores))
  results <- list()
  for (chunk_i in seq_along(chunks)) {
    genera <- chunks[[chunk_i]]
    idx <- match(genera, sort(unique(terminal$genus)))
    chunk_results <- parallel::mclapply(
      seq_along(genera),
      function(k) fit_fold(spec, cache, genera[[k]], idx[[k]]),
      mc.cores = min(cores, length(genera))
    )
    results <- c(results, chunk_results)
    fold_df <- do.call(rbind, lapply(results, `[[`, "fold"))
    checkpoint <- ifelse(smoke, paste0(spec$trait, ".smoke.fold_diagnostics.checkpoint.tsv"), paste0(spec$trait, ".fold_diagnostics.checkpoint.tsv"))
    write_tsv(fold_df, file.path(out_dir, checkpoint))
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
        species = excluded$Species,
        genus = excluded$genus,
        trait_category = excluded$trait_category,
        marine_binary = excluded$marine_binary,
        binary_aquatic_endpoint = NA_real_,
        aquaticity_score = excluded$aquaticity_score,
        held_out_genus = excluded$genus,
        globalFDR_corrected_OOF_prediction = NA_real_,
        nested_ttest_OOF_prediction = NA_real_,
        prediction_delta_nested_minus_globalFDR = NA_real_,
        nested_OOF_available = FALSE,
        exclusion_reason = "semi-aquatic aquatic_v2=0.5 excluded from binary aquatic-dependence LASSO/AUC",
        fold_n_ttest_genes_tested = NA_integer_,
        fold_n_candidate_genes_from_training_ttest = NA_integer_,
        fold_n_features_dropped_all_missing = NA_integer_,
        fold_n_features_dropped_zero_variance = NA_integer_,
        fold_n_features_used = NA_integer_,
        fold_n_selected_predictors = NA_integer_,
        fold_lambda = NA_real_,
        fold_model_status = "excluded_no_binary_aquatic_endpoint",
        stringsAsFactors = FALSE
      )
      oof <- rbind(oof, excluded_rows)
    }
    phase11_sub <- phase11_oof[phase11_oof$trait == spec$trait, c("species", "corrected_OOF_prediction"), drop = FALSE]
    names(phase11_sub)[2] <- "globalFDR_corrected_OOF_prediction"
    oof$globalFDR_corrected_OOF_prediction <- phase11_sub$globalFDR_corrected_OOF_prediction[match(oof$species, phase11_sub$species)]
    oof$prediction_delta_nested_minus_globalFDR <- oof$nested_ttest_OOF_prediction - oof$globalFDR_corrected_OOF_prediction
    oof <- oof[order(oof$species), , drop = FALSE]
  }
  eval_oof <- oof[oof$nested_OOF_available, , drop = FALSE]
  y <- if (spec$trait == "marine_binary") eval_oof$marine_binary else eval_oof$binary_aquatic_endpoint
  auc <- auc_rank(y, eval_oof$nested_ttest_OOF_prediction)
  list(spec = spec, fold_diag = fold_diag, oof = oof, selected = selected, warnings = warnings, auc = auc)
}

if (mode == "smoke") {
  stamp("Starting Phase 12B smoke test")
  smoke_results <- Map(function(spec) run_model(spec, trait_caches[[spec$trait]], smoke = TRUE), model_specs)
  smoke_folds <- do.call(rbind, lapply(smoke_results, `[[`, "fold_diag"))
  smoke_oof <- do.call(rbind, lapply(smoke_results, `[[`, "oof"))
  smoke_warnings <- do.call(rbind, lapply(smoke_results, `[[`, "warnings"))
  write_tsv(smoke_folds, file.path(smoke_log_dir, "nested_ttest_smoke_fold_details.tsv"))
  write_tsv(smoke_oof, file.path(smoke_log_dir, "nested_ttest_smoke_OOF_predictions.tsv"))
  if (is.null(smoke_warnings) || nrow(smoke_warnings) == 0) {
    smoke_warnings <- data.frame(trait = character(), model = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())
  }
  write_tsv(smoke_warnings, file.path(smoke_log_dir, "nested_ttest_smoke_warning_log.tsv"))
  smoke_prediction_summary <- aggregate(
    nested_ttest_OOF_prediction ~ model + held_out_genus,
    data = smoke_oof,
    FUN = function(x) paste(sprintf("%.4f", x), collapse = ",")
  )
  names(smoke_prediction_summary)[3] <- "held_out_predictions"
  smoke_summary <- merge(
    smoke_folds,
    smoke_prediction_summary,
    by = c("model", "held_out_genus"),
    all.x = TRUE,
    sort = FALSE
  )
  smoke_summary <- data.frame(
    model = smoke_summary$model,
    held_out_genus = smoke_summary$held_out_genus,
    n_training_terminal_species = smoke_summary$n_training_samples,
    n_test_terminal_species = smoke_summary$n_test_samples,
    n_training_positive = smoke_summary$n_train_positive,
    n_training_negative = smoke_summary$n_train_negative,
    n_candidate_features_FDR_0_01 = smoke_summary$n_fold_candidate_genes,
    n_features_after_foldwise_preprocessing = smoke_summary$n_features_used_for_lasso,
    n_selected_predictors = smoke_summary$n_selected_predictors,
    prediction_status = smoke_summary$fold_model_status,
    held_out_predictions = smoke_summary$held_out_predictions,
    runtime = sprintf("%.1f sec total smoke runtime", as.numeric(difftime(Sys.time(), script_start, units = "secs"))),
    status = ifelse(smoke_summary$n_test_samples > 0 & !is.na(smoke_summary$held_out_predictions), "PASS", "FAIL"),
    notes = smoke_summary$notes,
    stringsAsFactors = FALSE
  )
  write_tsv(smoke_summary, file.path(smoke_dir, "nested_ttest_smoke_test_summary.tsv"))
  writeLines(c(
    "# Phase 12B Nested-T-Test Smoke Test",
    "",
    "Smoke test completed before the full baseline run.",
    "",
    "Checks performed:",
    "- fold-specific training t-test/FDR candidate selection;",
    "- fold-wise training-only imputation and scaling;",
    "- terminal-only glmnet fit/lambda selection;",
    "- held-out genus prediction export.",
    "",
    "Smoke-test detailed fold logs are in `nested_ttest_smoke_test_logs/`.",
    "",
    paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
  ), file.path(smoke_dir, "nested_ttest_smoke_test_README.md"))
  stamp("Smoke test complete")
  quit(save = "no")
}

stamp("Starting full Phase 12B nested-t-test baseline robustness run with cores=", cores)
results <- Map(function(spec) run_model(spec, trait_caches[[spec$trait]], smoke = FALSE), model_specs)
fold_diagnostics <- do.call(rbind, lapply(results, `[[`, "fold_diag"))
all_oof <- do.call(rbind, lapply(results, `[[`, "oof"))
selected_by_fold <- do.call(rbind, lapply(results, `[[`, "selected"))
warning_log <- do.call(rbind, lapply(results, `[[`, "warnings"))
if (is.null(selected_by_fold)) selected_by_fold <- data.frame()
if (is.null(warning_log) || nrow(warning_log) == 0) {
  warning_log <- data.frame(trait = character(), model = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())
}

write_tsv(fold_diagnostics, file.path(out_dir, "nested_ttest_fold_candidate_gene_counts.tsv"))
oof_export <- data.frame(
  species = all_oof$species,
  genus = all_oof$genus,
  model = all_oof$model,
  trait_category = all_oof$trait_category,
  marine_binary = all_oof$marine_binary,
  binary_aquatic_endpoint = all_oof$binary_aquatic_endpoint,
  aquaticity_score = all_oof$aquaticity_score,
  held_out_genus = all_oof$held_out_genus,
  corrected_globalFDR_OOF_prediction = all_oof$globalFDR_corrected_OOF_prediction,
  nested_ttest_OOF_prediction = all_oof$nested_ttest_OOF_prediction,
  prediction_delta = all_oof$prediction_delta_nested_minus_globalFDR,
  nested_OOF_available = all_oof$nested_OOF_available,
  exclusion_reason = all_oof$exclusion_reason,
  fold_n_candidate_features_FDR = all_oof$fold_n_candidate_genes_from_training_ttest,
  fold_n_features_after_preprocessing = all_oof$fold_n_features_used,
  fold_n_selected_predictors = all_oof$fold_n_selected_predictors,
  fold_model_status = all_oof$fold_model_status,
  stringsAsFactors = FALSE
)
write_tsv(oof_export, file.path(out_dir, "nested_ttest_OOF_predictions.tsv"))
write_tsv(selected_by_fold, file.path(out_dir, "nested_ttest_fold_selected_predictors.tsv"))
write_tsv(warning_log, file.path(out_dir, "nested_ttest_warning_log.tsv"))
write_tsv(do.call(rbind, lapply(trait_caches, `[[`, "global_validation")), file.path(qc_dir, "nested_ttest_global_vectorized_validation.tsv"))

auc_summary <- do.call(rbind, lapply(results, function(res) {
  spec <- res$spec
  fold <- res$fold_diag
  oof_eval <- res$oof[res$oof$nested_OOF_available, , drop = FALSE]
  y <- if (spec$trait == "marine_binary") oof_eval$marine_binary else oof_eval$binary_aquatic_endpoint
  phase11_row <- phase11_auc[phase11_auc$trait == spec$trait, , drop = FALSE]
  data.frame(
    model = spec$model,
    trait_axis = spec$trait_axis,
    n_species_evaluated = sum(res$oof$nested_OOF_available),
    n_positive = sum(y == 1, na.rm = TRUE),
    n_negative = sum(y == 0, na.rm = TRUE),
    semi_aquatic_handling = spec$semi_aquatic_handling,
    AUC_legacy_global_preprocessing = phase11_row$AUC_legacy_global_preprocessing[[1]],
    AUC_corrected_globalFDR_foldwise_preprocessing = phase11_row$AUC_corrected_foldwise_preprocessing[[1]],
    AUC_nested_ttest = res$auc,
    delta_vs_corrected_globalFDR = res$auc - phase11_row$AUC_corrected_foldwise_preprocessing[[1]],
    n_folds = sum(fold$n_test_samples > 0),
    n_folds_zero_candidate_features = sum(fold$n_test_samples > 0 & fold$n_fold_candidate_genes == 0, na.rm = TRUE),
    n_folds_null_intercept_only = sum(fold$n_test_samples > 0 & grepl("^null_intercept_only", fold$fold_model_status)),
    median_candidate_features_per_fold = stats::median(fold$n_fold_candidate_genes[fold$n_test_samples > 0], na.rm = TRUE),
    IQR_candidate_features_per_fold = iqr_text(fold$n_fold_candidate_genes[fold$n_test_samples > 0], digits = 0),
    median_selected_predictors_per_fold = stats::median(fold$n_selected_predictors[fold$n_test_samples > 0], na.rm = TRUE),
    IQR_selected_predictors_per_fold = iqr_text(fold$n_selected_predictors[fold$n_test_samples > 0], digits = 0),
    evidence_level = "terminal-species genus-level LOOCV with supervised t-test/FDR candidate discovery repeated inside each outer training fold; global deterministic ASR branch states retained; fold-wise training-only imputation/scaling/lambda/model fitting",
    trait = spec$trait,
    comparison_model = phase11_row$model[[1]],
    stringsAsFactors = FALSE
  )
}))
write_tsv(auc_summary, file.path(out_dir, "nested_ttest_AUC_summary.tsv"))

comparison <- auc_summary[, c(
  "model", "trait_axis", "comparison_model", "AUC_legacy_global_preprocessing",
  "AUC_corrected_globalFDR_foldwise_preprocessing", "AUC_nested_ttest", "delta_vs_corrected_globalFDR",
  "n_species_evaluated", "median_candidate_features_per_fold", "IQR_candidate_features_per_fold", "evidence_level"
)]
write_tsv(comparison, file.path(out_dir, "nested_ttest_vs_globalFDR_AUC_comparison.tsv"))

category_summary <- do.call(rbind, lapply(split(all_oof, all_oof$trait), function(df) {
  do.call(rbind, lapply(split(df, df$trait_category), function(cat_df) {
    eval_df <- cat_df[cat_df$nested_OOF_available, , drop = FALSE]
    data.frame(
      trait = unique(df$trait),
      model = unique(df$model),
      trait_category = unique(cat_df$trait_category),
      n_total_rows = nrow(cat_df),
      n_nested_OOF_available = nrow(eval_df),
      globalFDR_corrected_median = ifelse(nrow(eval_df) > 0, stats::median(eval_df$globalFDR_corrected_OOF_prediction, na.rm = TRUE), NA_real_),
      nested_ttest_median = ifelse(nrow(eval_df) > 0, stats::median(eval_df$nested_ttest_OOF_prediction, na.rm = TRUE), NA_real_),
      nested_ttest_IQR = ifelse(nrow(eval_df) > 0, iqr_text(eval_df$nested_ttest_OOF_prediction, digits = 3), ""),
      nested_ttest_max = ifelse(nrow(eval_df) > 0, max(eval_df$nested_ttest_OOF_prediction, na.rm = TRUE), NA_real_),
      n_excluded = sum(!cat_df$nested_OOF_available),
      exclusion_reason = paste(unique(cat_df$exclusion_reason[cat_df$exclusion_reason != ""]), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(category_summary, file.path(out_dir, "nested_ttest_vs_corrected_globalFDR_category_summary.tsv"))
if (file.exists(paths$phase11_category)) {
  file.copy(paths$phase11_category, file.path(out_dir, "corrected_globalFDR_OOF_category_summary.tsv"), overwrite = TRUE)
}

roc_all <- do.call(rbind, lapply(split(all_oof[all_oof$nested_OOF_available, , drop = FALSE], all_oof$model[all_oof$nested_OOF_available]), function(df) {
  y <- if (unique(df$trait) == "marine_binary") df$marine_binary else df$binary_aquatic_endpoint
  auc <- auc_rank(y, df$nested_ttest_OOF_prediction)
  roc_points(y, df$nested_ttest_OOF_prediction, unique(df$model), auc)
}))
write_tsv(roc_all, file.path(out_dir, "nested_ttest_ROC_points.tsv"))

save_plot_all <- function(plot, stem, width, height) {
  grDevices::pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  grDevices::png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
  grDevices::svg(file.path(fig_dir, paste0(stem, ".svg")), width = width, height = height)
  print(plot)
  grDevices::dev.off()
}

roc_plot_df <- merge(
  roc_all,
  auc_summary[, c("model", "AUC_nested_ttest")],
  by = "model",
  all.x = TRUE
)
roc_plot_df$label <- ifelse(
  grepl("marine_binary", roc_plot_df$model),
  paste0("Marine model (AUC = ", sprintf("%.3f", roc_plot_df$AUC_nested_ttest), ")"),
  paste0("Binary aquatic-dependence model (AUC = ", sprintf("%.3f", roc_plot_df$AUC_nested_ttest), ")")
)
p_roc <- ggplot(roc_plot_df, aes(x = FPR, y = TPR, color = label)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey65", linewidth = 0.4) +
  geom_step(linewidth = 0.9) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_color_manual(values = c("#2F72C4", "#E28416")) +
  labs(x = "False positive rate", y = "True positive rate", color = NULL, title = "Nested t-test baseline ROC") +
  theme_classic(base_size = 9) +
  theme(legend.position = c(0.62, 0.18), legend.background = element_rect(fill = "white", color = "grey85"), plot.title = element_text(face = "bold", size = 10))
save_plot_all(p_roc, "Figure2B_nested_ttest_ROC_candidate", 4.2, 3.4)

dist_df <- all_oof[all_oof$nested_OOF_available, , drop = FALSE]
dist_df$panel <- ifelse(
  dist_df$trait == "marine_binary",
  "Marine-model OOF probability",
  "Binary aquatic-dependence-model OOF probability"
)
dist_df$panel <- factor(dist_df$panel, levels = c("Marine-model OOF probability", "Binary aquatic-dependence-model OOF probability"))
dist_df$category_plot <- factor(dist_df$trait_category, levels = c("terrestrial", "semi-aquatic", "non-marine aquatic", "marine"))
p_dist <- ggplot(dist_df, aes(x = category_plot, y = nested_ttest_OOF_prediction, color = trait_category)) +
  geom_boxplot(outlier.shape = NA, width = 0.56, color = "grey35", fill = "white", linewidth = 0.35) +
  geom_jitter(width = 0.16, height = 0, size = 1.2, alpha = 0.75) +
  facet_wrap(~ panel, nrow = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25), expand = c(0.02, 0.02)) +
  scale_color_manual(values = c("terrestrial" = "#8C8C8C", "semi-aquatic" = "#D9B44A", "non-marine aquatic" = "#E28416", "marine" = "#2F72C4")) +
  labs(x = NULL, y = "Nested-t-test OOF probability", title = "Out-of-fold prediction distributions") +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 9),
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold", size = 10)
  ) +
  geom_text(
    data = data.frame(
      panel = factor("Binary aquatic-dependence-model OOF probability", levels = c("Marine-model OOF probability", "Binary aquatic-dependence-model OOF probability")),
      x = 2.0, y = 0.085,
      label = "Semi-aquatic excluded\nfrom binary model/AUC"
    ),
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    size = 2.4,
    color = "grey30"
  )
save_plot_all(p_dist, "Figure2C_nested_ttest_OOF_distribution_candidate", 7.2, 3.4)

p_combined <- gridExtra::arrangeGrob(p_roc, p_dist, ncol = 2, widths = c(0.82, 1.35))
grDevices::pdf(file.path(fig_dir, "Figure2_BC_nested_ttest_candidate_combined.pdf"), width = 10.5, height = 3.7, useDingbats = FALSE)
grid::grid.draw(p_combined)
grDevices::dev.off()
grDevices::png(file.path(fig_dir, "Figure2_BC_nested_ttest_candidate_combined.png"), width = 10.5, height = 3.7, units = "in", res = 300)
grid::grid.draw(p_combined)
grDevices::dev.off()
grDevices::svg(file.path(fig_dir, "Figure2_BC_nested_ttest_candidate_combined.svg"), width = 10.5, height = 3.7)
grid::grid.draw(p_combined)
grDevices::dev.off()

design_audit <- data.frame(
  check = c(
    "GBI_rerun",
    "ASR_rerun",
    "nested_training_ttest_FDR_feature_discovery",
    "global_deterministic_ASR_branch_states_retained",
    "terminal_only_LASSO_samples",
    "held_out_terminal_branches_removed_from_fold_ttest",
    "foldwise_training_only_imputation",
    "foldwise_training_only_scaling",
    "training_only_lambda_selection",
    "glmnet_standardize_FALSE",
    "semi_aquatic_used_in_binary_aquatic_model"
  ),
  status = c(
    "PASS_no",
    "PASS_no",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_yes",
    "PASS_no"
  ),
  notes = c(
    "Endpoint-fix old-label no-fuse GBI matrix was reused.",
    "Existing deterministic ASR branch-state table was reused; no ASR rerun was performed.",
    "Candidate genes are selected separately within each outer held-out genus fold using training branch data.",
    "This robustness run isolates supervised feature-screening risk while retaining the accepted deterministic ASR layer.",
    "glmnet receives terminal species rows only.",
    "Held-out genus terminal branches are subtracted from the fold t-test/FDR discovery statistics.",
    "Missing values are imputed from training terminal samples only.",
    "Scaling parameters are estimated from imputed training terminal samples only.",
    "cv.glmnet uses training terminal samples only.",
    "Manual scaling is used before glmnet/cv.glmnet.",
    "aquatic_v2=0.5 species are excluded from binary aquatic training/testing/preprocessing/AUC."
  ),
  stringsAsFactors = FALSE
)
write_tsv(design_audit, file.path(qc_dir, "nested_ttest_design_audit.tsv"))

boundary_qc <- data.frame(
  check = c(
    "held_out_genus_used_for_ttest_feature_discovery",
    "held_out_genus_used_for_imputation",
    "held_out_genus_used_for_scaling",
    "held_out_genus_used_for_lambda",
    "held_out_genus_used_for_glmnet_fit",
    "internal_branches_used_for_lasso",
    "internal_branches_used_for_imputation_or_scaling",
    "internal_branches_used_for_training_ttest"
  ),
  status = c("PASS_no", "PASS_no", "PASS_no", "PASS_no", "PASS_no", "PASS_no", "PASS_no", "PASS_by_design_yes"),
  notes = c(
    "Held-out terminal genus branches are removed from fold t-test/FDR discovery; internal ASR branches are retained by design.",
    "Fold imputation means are computed from training terminal samples only.",
    "Fold scaling means/sds are computed from training terminal samples only.",
    "cv.glmnet receives training terminal rows only.",
    "glmnet receives training terminal rows only.",
    "No internal branches are included as LASSO samples.",
    "No internal branches are included in LASSO preprocessing.",
    "Internal branches contribute to branch-level discovery because the accepted design retains global deterministic ASR branch-level discovery."
  ),
  stringsAsFactors = FALSE
)
write_tsv(boundary_qc, file.path(qc_dir, "nested_ttest_boundary_check.tsv"))
scope_check <- data.frame(
  check = c(
    "GBI_ASR_frozen",
    "feature_selection_nested_inside_outer_fold",
    "held_out_genus_terminal_branches_excluded_from_ttest",
    "same_BH_FDR_rule_as_Stage1",
    "same_t_direction_as_Stage1",
    "foldwise_imputation_training_only",
    "foldwise_scaling_training_only",
    "glmnet_standardize_FALSE",
    "held_out_genus_not_used_for_lambda",
    "held_out_genus_not_used_for_glmnet_fit",
    "semi_aquatic_excluded_from_binary_aquatic_model",
    "drop_sensitivity_not_run",
    "official_Fig2_not_overwritten"
  ),
  status = c("yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes", "yes"),
  notes = c(
    "GBI and deterministic ASR branch-state table were reused; no GBI/ASR rerun.",
    "t-test/FDR candidate genes selected separately inside each outer held-out-genus fold.",
    "held-out genus terminal branch columns are subtracted from fold-level t-test statistics.",
    "old Stage 1 while-loop BH/FDR rule retained.",
    "tvalue = mean(state0) - mean(state1), matching Stage 1.",
    "imputation means are computed from training terminal samples only.",
    "scaling means/sds are computed from training terminal samples only.",
    "manual scaling is used and glmnet/cv.glmnet standardize=FALSE.",
    "cv.glmnet receives training terminal rows only.",
    "glmnet receives training terminal rows only.",
    "aquatic_v2=0.5 rows excluded from binary aquatic model/AUC.",
    "Phase 12B baseline only; no drop/sensitivity LASSO.",
    "Only provisional candidate Fig. 2B/C files were written in the Phase 12B folder."
  ),
  stringsAsFactors = FALSE
)
write_tsv(scope_check, file.path(qc_dir, "nested_ttest_scope_check.tsv"))
write_tsv(scope_check[scope_check$check %in% c("same_BH_FDR_rule_as_Stage1", "same_t_direction_as_Stage1"), , drop = FALSE], file.path(qc_dir, "nested_ttest_FDR_rule_check.tsv"))
write_tsv(boundary_qc, file.path(qc_dir, "nested_ttest_preprocessing_boundary_check.tsv"))

value_check <- data.frame(
  metric = c(
    "marine_AUC_globalFDR_corrected_preprocessing",
    "marine_AUC_nested_ttest",
    "marine_delta_AUC_nested_minus_globalFDR",
    "aquatic_AUC_globalFDR_corrected_preprocessing",
    "aquatic_AUC_nested_ttest",
    "aquatic_delta_AUC_nested_minus_globalFDR",
    "marine_species_evaluated",
    "aquatic_species_evaluated",
    "aquatic_semi_aquatic_excluded_rows"
  ),
  value = c(
    auc_summary$AUC_corrected_globalFDR_foldwise_preprocessing[auc_summary$trait == "marine_binary"],
    auc_summary$AUC_nested_ttest[auc_summary$trait == "marine_binary"],
    auc_summary$delta_vs_corrected_globalFDR[auc_summary$trait == "marine_binary"],
    auc_summary$AUC_corrected_globalFDR_foldwise_preprocessing[auc_summary$trait == "binary_aquatic_dependence"],
    auc_summary$AUC_nested_ttest[auc_summary$trait == "binary_aquatic_dependence"],
    auc_summary$delta_vs_corrected_globalFDR[auc_summary$trait == "binary_aquatic_dependence"],
    auc_summary$n_species_evaluated[auc_summary$trait == "marine_binary"],
    auc_summary$n_species_evaluated[auc_summary$trait == "binary_aquatic_dependence"],
    sum(all_oof$trait == "binary_aquatic_dependence" & !all_oof$nested_OOF_available)
  ),
  stringsAsFactors = FALSE
)
write_tsv(value_check, file.path(qc_dir, "nested_ttest_value_check.tsv"))

classify <- function(delta, auc) {
  if (!is.finite(delta) || !is.finite(auc)) return("STOP_EXECUTION_ERROR")
  if (abs(delta) <= 0.02 && auc >= 0.80) return("NESTED_PASS_STRONG")
  if (abs(delta) <= 0.05 && auc >= 0.70) return("NESTED_PASS_WITH_DOWNGRADE")
  "NESTED_WEAKENS_LASSO_MAIN"
}
decision_rows <- data.frame(
  trait = auc_summary$trait,
  AUC_corrected_globalFDR_foldwise_preprocessing = auc_summary$AUC_corrected_globalFDR_foldwise_preprocessing,
  AUC_nested_ttest = auc_summary$AUC_nested_ttest,
  delta_vs_corrected_globalFDR = auc_summary$delta_vs_corrected_globalFDR,
  decision_category = mapply(classify, auc_summary$delta_vs_corrected_globalFDR, auc_summary$AUC_nested_ttest),
  stringsAsFactors = FALSE
)
write_tsv(decision_rows, file.path(decision_dir, "nested_ttest_Result2.2_decision_categories.tsv"))

overall_decision <- if (any(decision_rows$decision_category == "STOP_EXECUTION_ERROR")) {
  "STOP_EXECUTION_ERROR"
} else if (any(decision_rows$decision_category == "NESTED_WEAKENS_LASSO_MAIN")) {
  "NESTED_WEAKENS_LASSO_MAIN"
} else if (any(decision_rows$decision_category == "NESTED_PASS_WITH_DOWNGRADE")) {
  "NESTED_PASS_WITH_DOWNGRADE"
} else {
  "NESTED_PASS_STRONG"
}

memo <- c(
  "# Phase 12B Nested-T-Test Baseline Robustness Decision Memo",
  "",
  paste0("Overall category: `", overall_decision, "`."),
  "",
  "## Design",
  "",
  "This robustness run tests whether supervised global branch-level t-test/FDR feature discovery materially inflates the baseline LASSO/gLOOCV evidence.",
  "",
  "Frozen aspects:",
  "- endpoint-fix old-label no-fuse GBI matrix retained;",
  "- existing deterministic ASR branch-state table retained;",
  "- terminal-only LASSO retained;",
  "- fold-wise training-only imputation/scaling/lambda/model fitting retained.",
  "",
  "Changed aspect:",
  "- t-test/FDR candidate-gene discovery is repeated inside each outer held-out-genus fold after removing that genus's terminal branches from the branch-level t-test statistics.",
  "",
  "This is a supervised feature-selection robustness audit, not a new GBI/ASR/enrichment run and not a drop/sensitivity rerun.",
  "",
  "## AUC comparison",
  "",
  paste(capture.output(print(decision_rows, row.names = FALSE)), collapse = "\n"),
  "",
  "## Interpretation",
  "",
  if (overall_decision == "NESTED_PASS_STRONG") {
    "The nested supervised t-test robustness run gives baseline AUCs close to the Phase 11 corrected global-FDR preprocessing results. Result 2.2/Fig. 2B-C can remain in the main evidence chain if described at the correct evidence level."
  } else if (overall_decision == "NESTED_PASS_WITH_DOWNGRADE") {
    "The nested supervised t-test robustness run preserves the baseline signal but changes magnitude enough that Result 2.2/Fig. 2B-C should use conservative wording and updated values."
  } else if (overall_decision == "NESTED_WEAKENS_LASSO_MAIN") {
    "The nested supervised t-test robustness run materially weakens one or both baseline LASSO results. Fig. 2B-C and Result 2.2 require major reframing before manuscript integration."
  } else {
    "Execution did not complete cleanly; do not interpret baseline LASSO robustness until errors are resolved."
  },
  "",
  "## Required decision questions",
  "",
  paste0("- Can Fig. 2B remain main? ", ifelse(overall_decision %in% c("NESTED_PASS_STRONG", "NESTED_PASS_WITH_DOWNGRADE"), "Yes, with nested-t-test values and explicit evidence-level wording.", "No; it requires major reframing before main-figure use.")),
  paste0("- Can Fig. 2C remain main? ", ifelse(overall_decision %in% c("NESTED_PASS_STRONG", "NESTED_PASS_WITH_DOWNGRADE"), "Yes, if redrawn from nested-t-test OOF predictions and labelled as model-behaviour visualization.", "No; current Fig. 2C evidence should be demoted or redesigned.")),
  paste0("- Can Result 2.2 retain 'predictive genomic profiles'? ", ifelse(overall_decision == "NESTED_PASS_STRONG", "Yes, but not as end-to-end feature-discovery prediction; use nested supervised t-test robustness wording.", ifelse(overall_decision == "NESTED_PASS_WITH_DOWNGRADE", "Only with conservative wording.", "No."))),
  paste0("- Does binary aquatic-dependence prediction survive? ", ifelse(decision_rows$decision_category[decision_rows$trait == "binary_aquatic_dependence"] %in% c("NESTED_PASS_STRONG", "NESTED_PASS_WITH_DOWNGRADE"), "Yes under the nested-t-test baseline robustness criteria.", "No or not without major downgrade.")),
  "- Should corrected-globalFDR AUCs be replaced, demoted, or shown only as conditional reference? Use nested-t-test AUCs for robustness-facing Fig. 2 decisions; keep corrected-globalFDR values only as conditional reference/comparison unless control chat decides otherwise.",
  "- Should Fig. 4/Fig. 5 LASSO architecture proceed under corrected-globalFDR design, or wait for nested-t-test sensitivity? This Phase 12B run is baseline only. Fig. 4/Fig. 5 architecture should not use legacy LASSO numbers; whether to extend nested t-test into sensitivity runs is a control-chat/Pro decision after reviewing this baseline result.",
  "",
  "## Phase 12A status",
  "",
  "Phase 12A corrected-preprocessing sensitivity/architecture work is deferred by user priority correction. Any scaffold created before this decision is provisional and was not mixed into this Phase 12B package.",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(memo, file.path(out_dir, "nested_ttest_Result2.2_decision_memo.md"))

run_summary <- c(
  "# Phase 12B Run Summary",
  "",
  paste0("Mode: full"),
  paste0("Cores requested: ", cores),
  paste0("Alpha/FDR threshold: ", alpha),
  "",
  "Outputs generated:",
  "- `nested_ttest_AUC_summary.tsv`",
  "- `nested_ttest_OOF_predictions.tsv`",
  "- `nested_ttest_fold_candidate_gene_counts.tsv`",
  "- `nested_ttest_vs_globalFDR_AUC_comparison.tsv`",
  "- `nested_ttest_Result2.2_decision_memo.md`",
  "",
  "No GBI, ASR, enrichment, manuscript text, final figures, or LASSO drop/sensitivity analyses were run.",
  "",
  "The run retained global deterministic ASR branch states and repeated supervised t-test/FDR candidate selection inside each outer held-out-genus fold by removing held-out terminal branches from fold-level t-test statistics.",
  "",
  paste0("Overall decision: `", overall_decision, "`."),
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(run_summary, file.path(qc_dir, "run_summary.md"))

readme <- c(
  "# Phase 12B Nested-T-Test Baseline Robustness",
  "",
  "This folder contains the baseline nested supervised t-test robustness audit requested after Phase 11.",
  "",
  "Purpose: quantify whether supervised global branch-level t-test/FDR feature discovery materially inflates Fig. 2B/C LASSO validation evidence.",
  "",
  "Baseline models run:",
  "- marine binary;",
  "- binary aquatic dependence.",
  "",
  "Frozen aspects:",
  "- endpoint-fix GBI retained;",
  "- deterministic ASR branch-state layer retained;",
  "- terminal-only LASSO retained;",
  "- fold-wise training-only imputation/scaling/lambda/model fitting retained.",
  "",
  "Changed aspect:",
  "- t-test/FDR candidate discovery repeated within each outer held-out-genus fold.",
  "",
  "Not run:",
  "- GBI;",
  "- ASR;",
  "- enrichment;",
  "- drop/sensitivity LASSO;",
  "- manuscript edits;",
  "- final figure generation.",
  "",
  paste0("Timestamp: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z"))
)
writeLines(readme, file.path(work_dir, "README_nested_ttest_baseline_robustness.md"))

stamp("Full Phase 12B nested-t-test baseline robustness run complete")
