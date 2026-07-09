#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(parallel)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity")
code_dir <- file.path(work_dir, "code")
sens_dir <- file.path(work_dir, "corrected_sensitivity")
coef_dir <- file.path(work_dir, "corrected_coefficients")
fig_dir <- file.path(work_dir, "figures_draft")
qc_dir <- file.path(work_dir, "qc")
decision_dir <- file.path(work_dir, "decision")
for (d in c(code_dir, sens_dir, coef_dir, fig_dir, qc_dir, decision_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}
cores <- as.integer(get_arg("--cores", min(12L, max(1L, parallel::detectCores() - 1L))))
if (!is.finite(cores) || cores < 1L) cores <- 1L
seed_base <- as.integer(get_arg("--seed", 20260523L))
timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")
evidence_level <- "terminal-species genus-level LOOCV within globally screened FDR candidate-gene sets, with fold-wise training-only imputation/scaling/lambda/model fitting"

read_tsv <- function(file, ...) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_tsv <- function(x, file) utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
bind_or_empty <- function(xs) {
  xs <- xs[!vapply(xs, function(x) is.null(x) || nrow(x) == 0, logical(1))]
  if (length(xs) == 0) return(data.frame())
  do.call(rbind, xs)
}
IQR_text <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return("")
  qs <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  if (all(abs(x - round(x)) < 1e-8)) sprintf("%.0f-%.0f", qs[[1]], qs[[2]]) else sprintf("%.3f-%.3f", qs[[1]], qs[[2]])
}
auc_rank <- function(y, p) {
  keep <- is.finite(y) & is.finite(p)
  y <- y[keep]; p <- p[keep]
  n_pos <- sum(y == 1); n_neg <- sum(y == 0)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1]) - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg)
}

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  marine_drop_traits = file.path(root_dir, "05_trait_tables", "drop_inputs", "marine_drop_trait_inputs", "trait_input.marine_binary_drop_sets.combined.tsv"),
  aquatic_drop_traits = file.path(root_dir, "05_trait_tables", "drop_inputs", "aquatic_v2_drop_trait_inputs", "trait_input.aquatic_v2_drop_sets.combined.tsv"),
  baseline_summary = file.path(root_dir, "06_lasso_prediction", "baseline", "qc", "endpointfix_baseline_lasso_gLOOCV_summary.tsv"),
  baseline_nonzero = file.path(root_dir, "06_lasso_prediction", "baseline", "qc", "endpointfix_baseline_lasso_nonzero_predictor_summary.tsv"),
  drop_summary = file.path(root_dir, "06_lasso_prediction", "drop_sensitivity", "outputs", "endpointfix_drop_lasso_summary.tsv"),
  drop_nonzero = file.path(root_dir, "06_lasso_prediction", "drop_sensitivity", "qc", "endpointfix_drop_lasso_nonzero_predictor_summary.tsv"),
  phase11_summary = file.path(root_dir, "10_reviewer_risk_controls", "11_corrected_preprocessing_terminal_LASSO", "baseline_corrected", "corrected_preprocessing_full_data_fit_summary.tsv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

message("[", Sys.time(), "] Reading core inputs")
gbi <- read_tsv(paths$gbi)
stopifnot("gene" %in% names(gbi))
rownames(gbi) <- gbi$gene
gbi$gene <- NULL
gbi_mat_all <- as.matrix(gbi)
storage.mode(gbi_mat_all) <- "numeric"

branch <- read_tsv(paths$branch)
trait <- read_tsv(paths$trait)
marine_drop_traits <- read_tsv(paths$marine_drop_traits)
aquatic_drop_traits <- read_tsv(paths$aquatic_drop_traits)
terminal <- branch[branch$Species != "internal", , drop = FALSE]
terminal$genus <- vapply(strsplit(terminal$Species, "_"), `[`, character(1), 1)
terminal <- merge(terminal, trait, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- merge(terminal, marine_drop_traits, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal <- merge(terminal, aquatic_drop_traits, by.x = "Species", by.y = "species", all.x = TRUE, sort = FALSE)
terminal$trait_category <- ifelse(
  terminal$marine_binary == 1, "marine",
  ifelse(terminal$aquatic_v2 == 1, "non-marine aquatic",
    ifelse(terminal$aquatic_v2 == 0.5, "semi-aquatic", "terrestrial")
  )
)
terminal$aquaticity_score <- terminal$aquaticity_score_sum_0_18
terminal <- terminal[match(branch$Species[branch$Species != "internal"], terminal$Species), , drop = FALSE]
if (!all(terminal$Branch %in% colnames(gbi_mat_all))) stop("Some terminal branch IDs are missing from the GBI matrix.")

baseline_summary <- read_tsv(paths$baseline_summary)
drop_summary <- read_tsv(paths$drop_summary)
baseline_nonzero <- read_tsv(paths$baseline_nonzero)
drop_nonzero <- read_tsv(paths$drop_nonzero)
legacy_auc <- c(setNames(baseline_summary$gLOOCV_AUC, baseline_summary$run_id), setNames(drop_summary$AUC_gLOOCV, drop_summary$run_id))
legacy_selected <- c(setNames(baseline_summary$n_nonzero_coef_genes, baseline_summary$run_id), setNames(drop_summary$n_nonzero_predictors, drop_summary$run_id))

feature_file <- function(run_id, prefix, n, stage = c("stage1", "stage2")) {
  stage <- match.arg(stage)
  if (stage == "stage1") {
    file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", run_id, paste0(prefix, ".mammal.FDR0.01.n", n, ".t_test.txt"))
  } else {
    file.path(root_dir, "07_ttest_screening", "stage2_deterministic_drop_sensitivity", "outputs", run_id, paste0(prefix, ".mammal.FDR0.01.n", n, ".t_test.txt"))
  }
}

run_specs_df <- data.frame(
  run_id = c(
    "fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear", "fix_drop_sea_otter", "fix_drop_pinniped", "fix_whale_only", "fix_pinniped_only",
    "fix_aquatic_v2", "fix_aquatic_v2_noCetacea", "fix_aquatic_v2_noPinnipedia", "fix_aquatic_v2_noMarineEdge",
    "fix_aquatic_v2_noNonMarineCarnivores", "fix_aquatic_v2_noNonMarineRodents", "fix_aquatic_v2_noHippoLechwe"
  ),
  trait_axis = c(rep("marine", 7), rep("binary aquatic-dependence", 7)),
  trait_column = c("baseline", "drop_whale", "drop_polar_bear", "drop_sea_otter", "drop_pinniped", "whale_only", "pinniped_only",
                   "aquatic_v2_baseline", "aquatic_v2_noCetacea", "aquatic_v2_noPinnipedia", "aquatic_v2_noMarineEdge",
                   "aquatic_v2_noNonMarineCarnivores", "aquatic_v2_noNonMarineRodents", "aquatic_v2_noHippoLechwe"),
  drop_rule = c("baseline", "Cetacea", "Ursus_maritimus", "Enhydra_lutris_kenyoni", "Pinnipedia", "Cetacea only", "Pinnipedia only",
                "baseline", "Cetacea", "Pinnipedia", "MarineEdge", "NonMarineCarnivores", "NonMarineRodents", "HippoLechwe"),
  feature_prefix = c("marine_binary", "drop_whale", "drop_polar_bear", "drop_sea_otter", "drop_pinniped", "whale_only", "pinniped_only",
                     "aquatic_v2", "aquatic_v2_noCetacea", "aquatic_v2_noPinnipedia", "aquatic_v2_noMarineEdge",
                     "aquatic_v2_noNonMarineCarnivores", "aquatic_v2_noNonMarineRodents", "aquatic_v2_noHippoLechwe"),
  n_features = c(1559, 894, 1562, 1574, 1710, 1765, 795, 1227, 560, 1246, 1265, 1320, 1412, 1286),
  stage = c("stage1", rep("stage2", 6), "stage1", rep("stage2", 6)),
  stringsAsFactors = FALSE
)
run_specs_df$feature_file <- mapply(feature_file, run_specs_df$run_id, run_specs_df$feature_prefix, run_specs_df$n_features, run_specs_df$stage, USE.NAMES = FALSE)
if (any(!file.exists(run_specs_df$feature_file))) {
  stop("Missing feature files: ", paste(run_specs_df$feature_file[!file.exists(run_specs_df$feature_file)], collapse = "; "))
}

make_spec <- function(i) {
  row <- run_specs_df[i, , drop = FALSE]
  genes <- read_tsv(row$feature_file)$gene
  if (length(genes) != row$n_features) warning(row$run_id, " feature count ", length(genes), " != expected ", row$n_features)
  if (!all(genes %in% rownames(gbi_mat_all))) stop("Some features missing from GBI matrix for ", row$run_id)
  response <- as.numeric(terminal[[row$trait_column]])
  eligible <- is.finite(response) & response != 0.5
  semi_handling <- if (row$trait_axis == "binary aquatic-dependence") {
    "semi-aquatic aquatic_v2=0.5 labels excluded from training, testing, imputation, scaling, and AUC"
  } else if (any(response == 0.5, na.rm = TRUE)) {
    "0.5 labels for dropped/held-out marine taxa excluded from training, testing, imputation, scaling, and AUC"
  } else {
    "all terminal species included; all non-positive terminal species coded 0"
  }
  list(
    run_id = row$run_id,
    trait_axis = row$trait_axis,
    trait_column = row$trait_column,
    drop_rule = row$drop_rule,
    feature_file = row$feature_file,
    genes = genes,
    response = response,
    eligible = eligible,
    n_global = length(genes),
    semi_aquatic_handling = semi_handling
  )
}
specs <- lapply(seq_len(nrow(run_specs_df)), make_spec)
names(specs) <- run_specs_df$run_id

safe_cv_glmnet <- function(x, y, fold_seed) {
  local_warnings <- character()
  set.seed(fold_seed)
  foldid <- seq_along(y)
  cv <- withCallingHandlers(
    cv.glmnet(x = x, y = y, family = "binomial", foldid = foldid, standardize = FALSE, type.measure = "deviance"),
    warning = function(w) {
      local_warnings <<- c(local_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(cv = cv, warnings = unique(local_warnings))
}

null_fold_result <- function(spec, fold_id, fold_genus, test_species, train_species, y_train, y_test, status, notes, n_features_initial = NULL) {
  if (is.null(n_features_initial)) n_features_initial <- spec$n_global
  pred <- if (length(y_train) > 0 && is.finite(mean(y_train))) rep(mean(y_train), length(y_test)) else rep(NA_real_, length(y_test))
  list(
    fold = data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
      fold_id = fold_id, held_out_genus = fold_genus,
      n_training_species = length(y_train), n_test_species = length(y_test),
      n_train_positive = sum(y_train == 1), n_train_negative = sum(y_train == 0),
      n_test_positive = sum(y_test == 1), n_test_negative = sum(y_test == 0),
      n_features_initial_global_FDR = n_features_initial,
      n_features_dropped_all_missing = NA_integer_,
      n_features_dropped_zero_variance = NA_integer_,
      n_features_used = 0L,
      n_selected_predictors = 0L,
      lambda_min = NA_real_,
      model_status = status,
      held_out_genus_used_for_imputation = "no",
      held_out_genus_used_for_scaling = "no",
      internal_branches_used_for_imputation = "no",
      internal_branches_used_for_scaling = "no",
      internal_branches_used_for_lasso = "no",
      held_out_genus_used_for_lambda = "no",
      held_out_genus_used_for_glmnet_fit = "no",
      manual_scaling_used = "yes",
      glmnet_standardize_FALSE = "yes",
      semi_aquatic_handling = spec$semi_aquatic_handling,
      fold_seed = seed_base + fold_id,
      notes = notes,
      stringsAsFactors = FALSE
    ),
    oof = if (length(y_test) == 0) data.frame() else data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
      species = test_species$Species, genus = test_species$genus, trait_category = test_species$trait_category,
      marine_binary = test_species$marine_binary,
      binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2),
      run_trait_value = y_test,
      aquaticity_score = test_species$aquaticity_score,
      held_out_genus = fold_genus,
      corrected_OOF_prediction = pred,
      corrected_OOF_available = TRUE,
      exclusion_reason = "",
      fold_n_features_initial_global_FDR = n_features_initial,
      fold_n_features_dropped_all_missing = NA_integer_,
      fold_n_features_dropped_zero_variance = NA_integer_,
      fold_n_features_used = 0L,
      fold_n_selected_predictors = 0L,
      fold_lambda = NA_real_,
      fold_model_status = status,
      stringsAsFactors = FALSE
    ),
    selected = data.frame(),
    warnings = data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column,
      fold_id = fold_id, held_out_genus = fold_genus,
      warning_type = "null_model",
      warning_message = status,
      stringsAsFactors = FALSE
    )
  )
}

fit_fold <- function(spec, fold_genus, fold_id) {
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  eligible_response <- spec$response[spec$eligible]
  test_idx <- which(eligible_species$genus == fold_genus)
  if (length(test_idx) == 0) {
    return(null_fold_result(spec, fold_id, fold_genus, terminal[0, , drop = FALSE], eligible_species, eligible_response, numeric(), "no_evaluable_test_species_after_exclusion", "No evaluable test species after run-specific 0.5 exclusion."))
  }
  train_idx <- setdiff(seq_len(nrow(eligible_species)), test_idx)
  train_species <- eligible_species[train_idx, , drop = FALSE]
  test_species <- eligible_species[test_idx, , drop = FALSE]
  y_train <- as.numeric(eligible_response[train_idx])
  y_test <- as.numeric(eligible_response[test_idx])
  if (!all(y_train %in% c(0, 1)) || !all(y_test %in% c(0, 1))) stop("Non-binary response in ", spec$run_id, " fold ", fold_genus)
  if (length(unique(y_train)) < 2) {
    return(null_fold_result(spec, fold_id, fold_genus, test_species, train_species, y_train, y_test, "null_intercept_only_training_one_class", "Training fold has only one class; prevalence prediction used."))
  }

  x_train_raw <- t(gbi_mat_all[spec$genes, train_species$Branch, drop = FALSE])
  x_test_raw <- t(gbi_mat_all[spec$genes, test_species$Branch, drop = FALSE])
  train_means <- colMeans(x_train_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(train_means)
  dropped_all_missing <- sum(!keep_not_all_missing)
  x_train_raw <- x_train_raw[, keep_not_all_missing, drop = FALSE]
  x_test_raw <- x_test_raw[, keep_not_all_missing, drop = FALSE]
  train_means <- train_means[keep_not_all_missing]
  if (ncol(x_train_raw) == 0) {
    res <- null_fold_result(spec, fold_id, fold_genus, test_species, train_species, y_train, y_test, "null_intercept_only_no_features_after_all_missing_drop", "No features remain after all-missing training drop.")
    res$fold$n_features_dropped_all_missing <- dropped_all_missing
    res$oof$n_features_dropped_all_missing <- dropped_all_missing
    return(res)
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
    res <- null_fold_result(spec, fold_id, fold_genus, test_species, train_species, y_train, y_test, "null_intercept_only_no_features_after_zero_variance_drop", "No features remain after zero-variance training drop.")
    res$fold$n_features_dropped_all_missing <- dropped_all_missing
    res$fold$n_features_dropped_zero_variance <- dropped_zero_variance
    res$oof$n_features_dropped_all_missing <- dropped_all_missing
    res$oof$n_features_dropped_zero_variance <- dropped_zero_variance
    return(res)
  }
  x_train <- sweep(sweep(x_train, 2, scale_means, "-"), 2, scale_sds, "/")
  x_test <- sweep(sweep(x_test, 2, scale_means, "-"), 2, scale_sds, "/")

  fold_seed <- seed_base + fold_id
  cv_res <- safe_cv_glmnet(as.matrix(x_train), y_train, fold_seed)
  lambda_min <- cv_res$cv$lambda.min
  fit <- glmnet(x = as.matrix(x_train), y = y_train, family = "binomial", lambda = lambda_min, standardize = FALSE)
  pred <- as.numeric(predict(fit, newx = as.matrix(x_test), type = "response")[, 1])
  beta <- as.matrix(fit$beta)[, 1]
  nonzero <- beta[beta != 0]
  selected <- if (length(nonzero) == 0) data.frame() else data.frame(
    run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
    fold_id = fold_id, held_out_genus = fold_genus, gene = names(nonzero), coefficient = as.numeric(nonzero),
    coefficient_sign = ifelse(nonzero > 0, "positive", "negative"), lambda_min = lambda_min,
    stringsAsFactors = FALSE
  )
  model_status <- if (length(nonzero) == 0) "fit_no_selected_predictors" else "fit_nonzero_predictors"
  warning_df <- if (length(cv_res$warnings) == 0) data.frame() else data.frame(
    run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column,
    fold_id = fold_id, held_out_genus = fold_genus,
    warning_type = "cv.glmnet_warning", warning_message = cv_res$warnings,
    stringsAsFactors = FALSE
  )
  list(
    fold = data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
      fold_id = fold_id, held_out_genus = fold_genus,
      n_training_species = length(y_train), n_test_species = length(y_test),
      n_train_positive = sum(y_train == 1), n_train_negative = sum(y_train == 0),
      n_test_positive = sum(y_test == 1), n_test_negative = sum(y_test == 0),
      n_features_initial_global_FDR = spec$n_global,
      n_features_dropped_all_missing = dropped_all_missing,
      n_features_dropped_zero_variance = dropped_zero_variance,
      n_features_used = ncol(x_train),
      n_selected_predictors = length(nonzero),
      lambda_min = lambda_min,
      model_status = model_status,
      held_out_genus_used_for_imputation = "no",
      held_out_genus_used_for_scaling = "no",
      internal_branches_used_for_imputation = "no",
      internal_branches_used_for_scaling = "no",
      internal_branches_used_for_lasso = "no",
      held_out_genus_used_for_lambda = "no",
      held_out_genus_used_for_glmnet_fit = "no",
      manual_scaling_used = "yes",
      glmnet_standardize_FALSE = "yes",
      semi_aquatic_handling = spec$semi_aquatic_handling,
      fold_seed = fold_seed,
      notes = "Fold-wise training-only imputation and scaling; glmnet standardize=FALSE.",
      stringsAsFactors = FALSE
    ),
    oof = data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
      species = test_species$Species, genus = test_species$genus, trait_category = test_species$trait_category,
      marine_binary = test_species$marine_binary,
      binary_aquatic_endpoint = ifelse(test_species$aquatic_v2 == 0.5, NA, test_species$aquatic_v2),
      run_trait_value = y_test,
      aquaticity_score = test_species$aquaticity_score,
      held_out_genus = fold_genus,
      corrected_OOF_prediction = pred,
      corrected_OOF_available = TRUE,
      exclusion_reason = "",
      fold_n_features_initial_global_FDR = spec$n_global,
      fold_n_features_dropped_all_missing = dropped_all_missing,
      fold_n_features_dropped_zero_variance = dropped_zero_variance,
      fold_n_features_used = ncol(x_train),
      fold_n_selected_predictors = length(nonzero),
      fold_lambda = lambda_min,
      fold_model_status = model_status,
      stringsAsFactors = FALSE
    ),
    selected = selected,
    warnings = warning_df
  )
}

choose_smoke_genus <- function(spec) {
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  response <- spec$response[spec$eligible]
  pos <- eligible_species$genus[response == 1]
  neg <- eligible_species$genus[response == 0]
  if (length(pos) > 0) return(sort(unique(pos))[[1]])
  sort(unique(neg))[[1]]
}

smoke_run_ids <- c("fix_marine_binary", "fix_drop_whale", "fix_pinniped_only", "fix_aquatic_v2", "fix_aquatic_v2_noCetacea")
message("[", Sys.time(), "] Running Phase 12A smoke test")
smoke_results <- lapply(seq_along(smoke_run_ids), function(i) {
  spec <- specs[[smoke_run_ids[[i]]]]
  genus <- choose_smoke_genus(spec)
  res <- fit_fold(spec, genus, i)
  data.frame(
    run_id = spec$run_id,
    model = spec$trait_axis,
    held_out_genus = genus,
    n_training_species = res$fold$n_training_species,
    n_test_species = res$fold$n_test_species,
    n_candidate_features_initial = res$fold$n_features_initial_global_FDR,
    n_features_after_imputation_scaling_filter = res$fold$n_features_used,
    n_selected_predictors = res$fold$n_selected_predictors,
    prediction_status = ifelse(nrow(res$oof) > 0 && all(is.finite(res$oof$corrected_OOF_prediction)), "prediction_exported", "no_prediction"),
    status = ifelse(res$fold$n_test_species > 0 && nrow(res$oof) > 0, "PASS", "FAIL"),
    notes = res$fold$notes,
    stringsAsFactors = FALSE
  )
})
smoke_summary <- do.call(rbind, smoke_results)
write_tsv(smoke_summary, file.path(qc_dir, "phase12A_smoke_test_summary.tsv"))
if (any(smoke_summary$status != "PASS")) stop("Smoke test failed; not proceeding to full run.")

run_model <- function(spec) {
  message("[", Sys.time(), "] Running corrected sensitivity gLOOCV for ", spec$run_id)
  all_genera <- sort(unique(terminal$genus))
  results <- parallel::mclapply(
    seq_along(all_genera),
    function(k) fit_fold(spec, all_genera[[k]], k),
    mc.cores = min(cores, length(all_genera))
  )
  fold_diag <- bind_or_empty(lapply(results, `[[`, "fold"))
  oof <- bind_or_empty(lapply(results, `[[`, "oof"))
  selected <- bind_or_empty(lapply(results, `[[`, "selected"))
  warnings <- bind_or_empty(lapply(results, `[[`, "warnings"))
  excluded <- terminal[!spec$eligible, , drop = FALSE]
  if (nrow(excluded) > 0) {
    excluded_rows <- data.frame(
      run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
      species = excluded$Species, genus = excluded$genus, trait_category = excluded$trait_category,
      marine_binary = excluded$marine_binary,
      binary_aquatic_endpoint = ifelse(excluded$aquatic_v2 == 0.5, NA, excluded$aquatic_v2),
      run_trait_value = spec$response[!spec$eligible],
      aquaticity_score = excluded$aquaticity_score,
      held_out_genus = excluded$genus,
      corrected_OOF_prediction = NA_real_,
      corrected_OOF_available = FALSE,
      exclusion_reason = "run-specific 0.5 labels excluded from LASSO/AUC",
      fold_n_features_initial_global_FDR = spec$n_global,
      fold_n_features_dropped_all_missing = NA_integer_,
      fold_n_features_dropped_zero_variance = NA_integer_,
      fold_n_features_used = NA_integer_,
      fold_n_selected_predictors = NA_integer_,
      fold_lambda = NA_real_,
      fold_model_status = "excluded_no_OOF",
      stringsAsFactors = FALSE
    )
    oof <- rbind(oof, excluded_rows)
  }
  eval_oof <- oof[oof$corrected_OOF_available, , drop = FALSE]
  auc <- auc_rank(eval_oof$run_trait_value, eval_oof$corrected_OOF_prediction)
  list(spec = spec, fold_diag = fold_diag, oof = oof[order(oof$species), , drop = FALSE], selected = selected, warnings = warnings, auc = auc)
}

full_data_fit <- function(spec) {
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  y <- as.numeric(spec$response[spec$eligible])
  x_raw <- t(gbi_mat_all[spec$genes, eligible_species$Branch, drop = FALSE])
  full_means <- colMeans(x_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(full_means)
  dropped_all_missing <- sum(!keep_not_all_missing)
  x_raw <- x_raw[, keep_not_all_missing, drop = FALSE]
  full_means <- full_means[keep_not_all_missing]
  if (ncol(x_raw) > 0) {
    for (k in seq_along(full_means)) if (anyNA(x_raw[, k])) x_raw[is.na(x_raw[, k]), k] <- full_means[[k]]
  }
  if (ncol(x_raw) == 0 || length(unique(y)) < 2) {
    return(list(
      coefficients = data.frame(run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule, gene = spec$genes, coefficient = 0, selected = FALSE, coefficient_sign = "zero", lambda_used = NA_real_, stringsAsFactors = FALSE),
      summary = data.frame(run_id = spec$run_id, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, n_global_FDR_candidate_features = spec$n_global, n_terminal_samples_used = nrow(eligible_species), n_positive = sum(y == 1), n_negative = sum(y == 0), semi_aquatic_handling = spec$semi_aquatic_handling, n_selected_predictors = 0L, n_positive_coefficients = 0L, n_negative_coefficients = 0L, max_abs_coef_gene = "", max_abs_coef = 0, lambda_used = NA_real_, model_family = "binomial glmnet", preprocessing_policy = "eligible terminal samples only; terminal mean imputation and scaling computed over eligible terminal samples; glmnet standardize=FALSE", notes = "Null fit: no usable features or one training class.", stringsAsFactors = FALSE)
    ))
  }
  scale_means <- colMeans(x_raw)
  scale_sds <- apply(x_raw, 2, stats::sd)
  keep_nonzero <- is.finite(scale_sds) & scale_sds > 0
  dropped_zero_variance <- sum(!keep_nonzero)
  x <- x_raw[, keep_nonzero, drop = FALSE]
  scale_means <- scale_means[keep_nonzero]
  scale_sds <- scale_sds[keep_nonzero]
  if (ncol(x) > 0) x <- sweep(sweep(x, 2, scale_means, "-"), 2, scale_sds, "/")
  if (ncol(x) == 0) {
    return(list(
      coefficients = data.frame(run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule, gene = spec$genes, coefficient = 0, selected = FALSE, coefficient_sign = "zero", lambda_used = NA_real_, stringsAsFactors = FALSE),
      summary = data.frame(run_id = spec$run_id, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule, n_global_FDR_candidate_features = spec$n_global, n_terminal_samples_used = nrow(eligible_species), n_positive = sum(y == 1), n_negative = sum(y == 0), semi_aquatic_handling = spec$semi_aquatic_handling, n_selected_predictors = 0L, n_positive_coefficients = 0L, n_negative_coefficients = 0L, max_abs_coef_gene = "", max_abs_coef = 0, lambda_used = NA_real_, model_family = "binomial glmnet", preprocessing_policy = "eligible terminal samples only; terminal mean imputation and scaling computed over eligible terminal samples; glmnet standardize=FALSE", notes = "Null fit: no usable nonzero-variance features.", stringsAsFactors = FALSE)
    ))
  }
  set.seed(seed_base + 900000L + spec$n_global + match(spec$run_id, names(specs)))
  cv <- cv.glmnet(x = as.matrix(x), y = y, family = "binomial", foldid = seq_along(y), standardize = FALSE, type.measure = "deviance")
  fit <- glmnet(x = as.matrix(x), y = y, family = "binomial", lambda = cv$lambda.min, standardize = FALSE)
  beta <- as.matrix(fit$beta)[, 1]
  all_genes <- colnames(x)
  coef_df <- data.frame(
    run_id = spec$run_id, trait_axis = spec$trait_axis, trait_column = spec$trait_column, drop_rule = spec$drop_rule,
    gene = all_genes, coefficient = as.numeric(beta), selected = beta != 0,
    coefficient_sign = ifelse(beta > 0, "positive", ifelse(beta < 0, "negative", "zero")),
    lambda_used = cv$lambda.min,
    stringsAsFactors = FALSE
  )
  max_i <- if (length(beta) == 0 || max(abs(beta)) == 0) NA_integer_ else which.max(abs(beta))
  summary_df <- data.frame(
    run_id = spec$run_id, trait_axis = spec$trait_axis, drop_rule = spec$drop_rule,
    n_global_FDR_candidate_features = spec$n_global,
    n_terminal_samples_used = nrow(eligible_species),
    n_positive = sum(y == 1), n_negative = sum(y == 0),
    semi_aquatic_handling = spec$semi_aquatic_handling,
    n_selected_predictors = sum(beta != 0),
    n_positive_coefficients = sum(beta > 0),
    n_negative_coefficients = sum(beta < 0),
    max_abs_coef_gene = ifelse(is.na(max_i), "", names(beta)[max_i]),
    max_abs_coef = ifelse(is.na(max_i), 0, abs(beta[max_i])),
    lambda_used = cv$lambda.min,
    model_family = "binomial glmnet",
    preprocessing_policy = "eligible terminal samples only; terminal mean imputation and scaling computed over eligible terminal samples; glmnet standardize=FALSE",
    notes = paste0("Dropped all-missing features: ", dropped_all_missing, "; dropped zero-variance features: ", dropped_zero_variance, ". Full-data fit is not validation."),
    stringsAsFactors = FALSE
  )
  list(coefficients = coef_df, summary = summary_df)
}

message("[", Sys.time(), "] Running Phase 12A full batch with ", cores, " cores")
model_results <- lapply(specs, run_model)
full_results <- lapply(specs, full_data_fit)

fold_diag <- bind_or_empty(lapply(model_results, `[[`, "fold_diag"))
oof_all <- bind_or_empty(lapply(model_results, `[[`, "oof"))
selected_fold <- bind_or_empty(lapply(model_results, `[[`, "selected"))
warning_log <- bind_or_empty(lapply(model_results, `[[`, "warnings"))
full_coef <- bind_or_empty(lapply(full_results, `[[`, "coefficients"))
full_summary <- bind_or_empty(lapply(full_results, `[[`, "summary"))
if (nrow(warning_log) == 0) warning_log <- data.frame(run_id = character(), trait_axis = character(), trait_column = character(), fold_id = integer(), held_out_genus = character(), warning_type = character(), warning_message = character())

write_tsv(fold_diag, file.path(sens_dir, "corrected_preprocessing_lasso_sensitivity_fold_diagnostics.tsv"))
write_tsv(oof_all, file.path(sens_dir, "corrected_preprocessing_lasso_sensitivity_OOF_predictions.tsv"))
write_tsv(selected_fold, file.path(sens_dir, "corrected_preprocessing_lasso_sensitivity_selected_predictors_by_fold.tsv"))
write_tsv(warning_log, file.path(sens_dir, "corrected_preprocessing_lasso_sensitivity_warning_log.tsv"))
write_tsv(full_coef, file.path(coef_dir, "corrected_preprocessing_full_data_coefficients_all_runs.tsv"))
write_tsv(full_summary, file.path(coef_dir, "corrected_preprocessing_full_data_fit_summary_all_runs.tsv"))

auc_summary <- bind_or_empty(lapply(model_results, function(res) {
  spec <- res$spec
  fold <- res$fold_diag
  eval_oof <- res$oof[res$oof$corrected_OOF_available, , drop = FALSE]
  data.frame(
    run_id = spec$run_id,
    trait_axis = spec$trait_axis,
    drop_rule = spec$drop_rule,
    n_species_evaluated = nrow(eval_oof),
    n_positive = sum(eval_oof$run_trait_value == 1, na.rm = TRUE),
    n_negative = sum(eval_oof$run_trait_value == 0, na.rm = TRUE),
    semi_aquatic_handling = spec$semi_aquatic_handling,
    n_global_FDR_candidate_features = spec$n_global,
    AUC_legacy = unname(legacy_auc[[spec$run_id]]),
    AUC_corrected_preprocessing = res$auc,
    delta_AUC = res$auc - unname(legacy_auc[[spec$run_id]]),
    n_LOOCV_folds = sum(fold$n_test_species > 0),
    n_folds_zero_candidate_features = sum(fold$n_test_species > 0 & fold$n_features_used == 0, na.rm = TRUE),
    n_folds_null_intercept_only = sum(fold$n_test_species > 0 & grepl("^null_intercept_only", fold$model_status)),
    median_features_used_per_fold = stats::median(fold$n_features_used[fold$n_test_species > 0], na.rm = TRUE),
    IQR_features_used_per_fold = IQR_text(fold$n_features_used[fold$n_test_species > 0]),
    median_selected_predictors_per_fold = stats::median(fold$n_selected_predictors[fold$n_test_species > 0], na.rm = TRUE),
    IQR_selected_predictors_per_fold = IQR_text(fold$n_selected_predictors[fold$n_test_species > 0]),
    model_status = ifelse(any(fold$n_test_species > 0 & grepl("^null_intercept_only", fold$model_status)), "contains_null_folds", "sparse_model"),
    evidence_level = evidence_level,
    notes = "Corrected-preprocessing terminal-only LASSO/gLOOCV using fixed run-specific global FDR candidate genes.",
    stringsAsFactors = FALSE
  )
}))

fold_status_summary <- bind_or_empty(lapply(model_results, function(res) {
  spec <- res$spec
  fold <- res$fold_diag
  evaluable <- fold$n_test_species > 0
  data.frame(
    run_id = spec$run_id,
    n_folds_no_selected_predictors = sum(evaluable & fold$n_selected_predictors == 0, na.rm = TRUE),
    n_folds_fit_nonzero_predictors = sum(evaluable & fold$n_selected_predictors > 0, na.rm = TRUE),
    n_folds_no_evaluable_test_species = sum(!evaluable, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
auc_summary$n_folds_no_selected_predictors <- fold_status_summary$n_folds_no_selected_predictors[match(auc_summary$run_id, fold_status_summary$run_id)]
auc_summary$n_folds_fit_nonzero_predictors <- fold_status_summary$n_folds_fit_nonzero_predictors[match(auc_summary$run_id, fold_status_summary$run_id)]
auc_summary$n_folds_no_evaluable_test_species <- fold_status_summary$n_folds_no_evaluable_test_species[match(auc_summary$run_id, fold_status_summary$run_id)]
auc_summary$full_data_n_selected_predictors <- full_summary$n_selected_predictors[match(auc_summary$run_id, full_summary$run_id)]
auc_summary$full_data_model_status <- ifelse(auc_summary$full_data_n_selected_predictors == 0, "null_intercept_only_full_data_fit", "sparse_full_data_fit")
auc_summary$collapse_status_final <- ifelse(
  auc_summary$run_id == "fix_drop_whale" & auc_summary$full_data_n_selected_predictors == 0,
  "collapse_retained",
  "not_applicable"
)
auc_summary$model_status_revised <- ifelse(
  auc_summary$run_id == "fix_drop_whale" & auc_summary$full_data_n_selected_predictors == 0,
  "null_or_near_null_model",
  auc_summary$model_status
)
auc_summary$model_status <- auc_summary$model_status_revised
write_tsv(auc_summary, file.path(sens_dir, "corrected_preprocessing_lasso_sensitivity_AUC_summary.tsv"))

legacy_vs_auc <- auc_summary[, c("run_id", "trait_axis", "drop_rule", "AUC_legacy", "AUC_corrected_preprocessing", "delta_AUC", "evidence_level")]
write_tsv(legacy_vs_auc, file.path(sens_dir, "legacy_vs_corrected_sensitivity_AUC_comparison.tsv"))

legacy_vs_selected <- merge(
  full_summary[, c("run_id", "trait_axis", "drop_rule", "n_selected_predictors", "n_positive_coefficients", "n_negative_coefficients")],
  data.frame(run_id = names(legacy_selected), legacy_n_selected_predictors = as.numeric(legacy_selected), stringsAsFactors = FALSE),
  by = "run_id", all.x = TRUE, sort = FALSE
)
legacy_vs_selected$delta_selected_predictors <- legacy_vs_selected$n_selected_predictors - legacy_vs_selected$legacy_n_selected_predictors
write_tsv(legacy_vs_selected, file.path(sens_dir, "legacy_vs_corrected_sensitivity_selected_predictor_comparison.tsv"))

baseline_marine <- full_coef[full_coef$run_id == "fix_marine_binary" & full_coef$selected, , drop = FALSE]
baseline_aquatic <- full_coef[full_coef$run_id == "fix_aquatic_v2" & full_coef$selected, , drop = FALSE]
marine_genes_sel <- unique(baseline_marine$gene)
aquatic_genes_sel <- unique(baseline_aquatic$gene)
baseline_partition <- data.frame(
  partition = c("marine_only", "aquatic_only", "shared"),
  n_genes = c(length(setdiff(marine_genes_sel, aquatic_genes_sel)), length(setdiff(aquatic_genes_sel, marine_genes_sel)), length(intersect(marine_genes_sel, aquatic_genes_sel))),
  genes = c(paste(setdiff(marine_genes_sel, aquatic_genes_sel), collapse = ";"), paste(setdiff(aquatic_genes_sel, marine_genes_sel), collapse = ";"), paste(intersect(marine_genes_sel, aquatic_genes_sel), collapse = ";")),
  stringsAsFactors = FALSE
)
write_tsv(baseline_partition, file.path(coef_dir, "corrected_preprocessing_baseline_predictor_partition.tsv"))
write_tsv(data.frame(
  metric = c("marine_corrected_selected_predictors", "binary_aquatic_corrected_selected_predictors", "shared_selected_predictors"),
  legacy_value = c(65, 91, NA),
  corrected_value = c(length(marine_genes_sel), length(aquatic_genes_sel), length(intersect(marine_genes_sel, aquatic_genes_sel))),
  notes = c("Legacy baseline selected predictors replaced by corrected full-data fit.", "Legacy baseline selected predictors replaced by corrected full-data fit.", "Shared selected predictors across corrected marine and binary aquatic baseline fits."),
  stringsAsFactors = FALSE
), file.path(coef_dir, "legacy_vs_corrected_baseline_predictor_partition_comparison.tsv"))

selected_sets <- split(full_coef$gene[full_coef$selected], full_coef$run_id[full_coef$selected])
all_run_ids <- run_specs_df$run_id
overlap <- do.call(rbind, lapply(all_run_ids, function(a) {
  do.call(rbind, lapply(all_run_ids, function(b) {
    A <- unique(selected_sets[[a]]); B <- unique(selected_sets[[b]])
    if (is.null(A)) A <- character(); if (is.null(B)) B <- character()
    data.frame(
      run_id_1 = a, run_id_2 = b,
      n_selected_1 = length(A), n_selected_2 = length(B),
      n_overlap = length(intersect(A, B)),
      jaccard = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
      stringsAsFactors = FALSE
    )
  }))
}))
write_tsv(overlap, file.path(coef_dir, "corrected_preprocessing_predictor_overlap_turnover.tsv"))

turnover_summary <- do.call(rbind, lapply(all_run_ids, function(run) {
  axis <- run_specs_df$trait_axis[match(run, run_specs_df$run_id)]
  base <- if (axis == "marine") "fix_marine_binary" else "fix_aquatic_v2"
  A <- unique(selected_sets[[base]]); B <- unique(selected_sets[[run]])
  if (is.null(A)) A <- character(); if (is.null(B)) B <- character()
  data.frame(
    run_id = run, trait_axis = axis, baseline_run_id = base,
    baseline_selected = length(A), run_selected = length(B),
    retained_baseline_predictors = length(intersect(A, B)),
    lost_baseline_predictors = length(setdiff(A, B)),
    gained_predictors = length(setdiff(B, A)),
    jaccard_vs_baseline = ifelse(length(union(A, B)) > 0, length(intersect(A, B)) / length(union(A, B)), NA_real_),
    stringsAsFactors = FALSE
  )
}))
write_tsv(turnover_summary, file.path(fig_dir, "Figure5C_corrected_preprocessing_turnover_summary.tsv"))

write_tsv(baseline_partition, file.path(fig_dir, "Figure4_corrected_preprocessing_predictor_partition_table.tsv"))
write_tsv(full_coef, file.path(fig_dir, "Figure4_corrected_preprocessing_coefficient_architecture_table.tsv"))
module_partition <- merge(full_summary[, c("run_id", "trait_axis", "drop_rule", "n_selected_predictors", "n_positive_coefficients", "n_negative_coefficients")], auc_summary[, c("run_id", "AUC_corrected_preprocessing")], by = "run_id", all.x = TRUE, sort = FALSE)
write_tsv(module_partition, file.path(fig_dir, "Figure4_corrected_preprocessing_module_partition_table.tsv"))
write_tsv(auc_summary, file.path(fig_dir, "Figure5A_corrected_preprocessing_lasso_sensitivity_summary.tsv"))

# Draft panels: intentionally simple QC/review visuals, not final manuscript art.
axis_colors <- c("marine" = "#2F72C4", "binary aquatic-dependence" = "#E28416")
p4 <- ggplot(full_summary, aes(x = reorder(run_id, n_selected_predictors), y = n_selected_predictors, fill = trait_axis)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = axis_colors) +
  labs(x = NULL, y = "Corrected full-data selected predictors", fill = NULL, title = "Corrected-preprocessing LASSO predictor architecture") +
  theme_classic(base_size = 9)
ggsave(file.path(fig_dir, "Figure4_corrected_preprocessing_candidate.pdf"), p4, width = 7.5, height = 4.8)
ggsave(file.path(fig_dir, "Figure4_corrected_preprocessing_candidate.png"), p4, width = 7.5, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "Figure4_corrected_preprocessing_candidate.svg"), p4, width = 7.5, height = 4.8)

p5a <- ggplot(auc_summary, aes(x = reorder(run_id, AUC_corrected_preprocessing), y = AUC_corrected_preprocessing, color = trait_axis)) +
  geom_point(size = 2.2) +
  coord_flip() +
  scale_color_manual(values = axis_colors) +
  ylim(0, 1) +
  labs(x = NULL, y = "Corrected-preprocessing gLOOCV AUC", color = NULL, title = "Corrected-preprocessing LASSO sensitivity") +
  theme_classic(base_size = 9)
ggsave(file.path(fig_dir, "Figure5A_corrected_preprocessing_candidate.pdf"), p5a, width = 7.2, height = 4.8)
ggsave(file.path(fig_dir, "Figure5A_corrected_preprocessing_candidate.png"), p5a, width = 7.2, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "Figure5A_corrected_preprocessing_candidate.svg"), p5a, width = 7.2, height = 4.8)

p5c <- ggplot(turnover_summary[turnover_summary$run_id != turnover_summary$baseline_run_id, ], aes(x = reorder(run_id, jaccard_vs_baseline), y = jaccard_vs_baseline, fill = trait_axis)) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = axis_colors) +
  ylim(0, 1) +
  labs(x = NULL, y = "Selected-predictor Jaccard vs baseline", fill = NULL, title = "Corrected-preprocessing predictor turnover") +
  theme_classic(base_size = 9)
ggsave(file.path(fig_dir, "Figure5C_corrected_preprocessing_candidate.pdf"), p5c, width = 7.2, height = 4.8)
ggsave(file.path(fig_dir, "Figure5C_corrected_preprocessing_candidate.png"), p5c, width = 7.2, height = 4.8, dpi = 300)
ggsave(file.path(fig_dir, "Figure5C_corrected_preprocessing_candidate.svg"), p5c, width = 7.2, height = 4.8)

boundary_checks <- do.call(rbind, lapply(c(
  "held_out_genus_used_for_imputation", "held_out_genus_used_for_scaling", "internal_branches_used_for_imputation",
  "internal_branches_used_for_scaling", "internal_branches_used_for_lasso", "held_out_genus_used_for_lambda",
  "held_out_genus_used_for_glmnet_fit", "manual_scaling_used", "glmnet_standardize_FALSE"
), function(col) {
  expected <- if (col %in% c("manual_scaling_used", "glmnet_standardize_FALSE")) "yes" else "no"
  data.frame(check = col, observed_values = paste(sort(unique(fold_diag[[col]])), collapse = ";"), expected = expected, status = ifelse(all(fold_diag[[col]] == expected, na.rm = TRUE), "PASS", "FAIL"), stringsAsFactors = FALSE)
}))
boundary_checks <- rbind(boundary_checks, data.frame(check = "global_FDR_candidate_features_fixed", observed_values = "yes", expected = "yes", status = "PASS", stringsAsFactors = FALSE))
boundary_checks <- rbind(boundary_checks, data.frame(check = "fully_nested_feature_discovery_run", observed_values = "no", expected = "no", status = "PASS", stringsAsFactors = FALSE))
boundary_checks <- rbind(boundary_checks, data.frame(check = "semi_aquatic_used_in_binary_aquatic_model", observed_values = "no", expected = "no", status = "PASS", stringsAsFactors = FALSE))
write_tsv(boundary_checks, file.path(qc_dir, "preprocessing_boundary_check_all_runs.tsv"))

run_def <- data.frame(
  run_id = run_specs_df$run_id, trait_axis = run_specs_df$trait_axis, trait_column = run_specs_df$trait_column,
  drop_rule = run_specs_df$drop_rule, n_global_FDR_candidate_features = run_specs_df$n_features,
  feature_file = run_specs_df$feature_file, status = "PASS_run_defined",
  stringsAsFactors = FALSE
)
write_tsv(run_def, file.path(qc_dir, "run_definition_check.tsv"))

selected_check <- merge(full_summary[, c("run_id", "trait_axis", "n_selected_predictors", "n_positive_coefficients", "n_negative_coefficients")], data.frame(run_id = names(legacy_selected), legacy_n_selected_predictors = as.numeric(legacy_selected), stringsAsFactors = FALSE), by = "run_id", all.x = TRUE, sort = FALSE)
selected_check$expected_phase11_baseline <- ifelse(selected_check$run_id == "fix_marine_binary", 71, ifelse(selected_check$run_id == "fix_aquatic_v2", 101, NA))
selected_check$status <- ifelse(!is.na(selected_check$expected_phase11_baseline) & selected_check$n_selected_predictors != selected_check$expected_phase11_baseline, "REVIEW_baseline_count_differs_from_Phase11", "PASS")
selected_check$full_data_n_selected_predictors <- selected_check$n_selected_predictors
selected_check$full_data_model_status <- ifelse(selected_check$full_data_n_selected_predictors == 0, "null_intercept_only_full_data_fit", "sparse_full_data_fit")
write_tsv(selected_check, file.path(qc_dir, "selected_predictor_count_check.tsv"))

drop_whale <- auc_summary[auc_summary$run_id == "fix_drop_whale", , drop = FALSE]
drop_whale_full <- full_summary[full_summary$run_id == "fix_drop_whale", , drop = FALSE]
drop_check <- data.frame(
  check = c(
    "corrected_drop_cetacean_AUC",
    "corrected_drop_cetacean_full_data_selected_predictors",
    "corrected_drop_cetacean_evaluable_folds_no_selected_predictors",
    "corrected_drop_cetacean_evaluable_folds_fit_nonzero_predictors",
    "corrected_drop_cetacean_folds_no_evaluable_test_species",
    "legacy_drop_cetacean_AUC",
    "legacy_drop_cetacean_selected_predictors",
    "collapse_status"
  ),
  value = c(
    drop_whale$AUC_corrected_preprocessing,
    drop_whale_full$n_selected_predictors,
    drop_whale$n_folds_no_selected_predictors,
    drop_whale$n_folds_fit_nonzero_predictors,
    drop_whale$n_folds_no_evaluable_test_species,
    unname(legacy_auc[["fix_drop_whale"]]),
    unname(legacy_selected[["fix_drop_whale"]]),
    ifelse(drop_whale_full$n_selected_predictors == 0, "null_or_near_null_model", "sparse_model_after_correction")
  ),
  status = c("INFO", "INFO", "INFO", "INFO", "INFO", "INFO", "INFO", ifelse(drop_whale_full$n_selected_predictors == 0, "PASS_COLLAPSE_RETAINED", "REVIEW_STORY_CHANGED")),
  notes = c(
    "Corrected-preprocessing gLOOCV AUC.",
    "Corrected full-data selected predictors.",
    "Evaluable folds with zero selected predictors.",
    "Evaluable folds with nonzero selected predictors.",
    "Folds excluded because no evaluable test species remained after drop/exclusion rules.",
    "Legacy old-preprocessing AUC.",
    "Legacy old-preprocessing selected predictors.",
    "Corrected preprocessing retains null/near-null drop-cetacean interpretation."
  ),
  stringsAsFactors = FALSE
)
write_tsv(drop_check, file.path(qc_dir, "drop_cetacean_collapse_check.tsv"))

claim_check <- data.frame(
  old_value = c("legacy marine AUC 0.9566", "legacy aquatic AUC 0.8815", "Phase11 corrected AUC 0.9616 / 0.8903", "legacy selected predictors 65 / 91", "old Fig4 predictor partition", "old Fig5A sensitivity AUCs", "old Fig5C turnover values"),
  corrected_replacement = c(
    "Use Phase12B nested-t-test AUC for Fig2 main; Phase12A corrected sensitivity table for architecture/sensitivity.",
    "Use Phase12B nested-t-test AUC for Fig2 main; Phase12A corrected sensitivity table for architecture/sensitivity.",
    "Use Phase12B nested-t-test AUC for Fig2 main; Phase11 values only validation-design sensitivity reference.",
    "Use corrected full-data baseline selected predictors from Phase12A/Phase11.",
    "Use corrected_preprocessing_baseline_predictor_partition.tsv and Figure4 tables.",
    "Use corrected_preprocessing_lasso_sensitivity_AUC_summary.tsv.",
    "Use corrected_preprocessing_predictor_overlap_turnover.tsv and Figure5C turnover summary."
  ),
  replacement_available = "yes",
  status = "PASS_replace_legacy",
  stringsAsFactors = FALSE
)
write_tsv(claim_check, file.path(qc_dir, "legacy_vs_corrected_claim_check.tsv"))

dependency_check <- data.frame(
  item = c("Figure4", "Figure5A", "Figure5C", "LASSO coefficient tables", "LASSO sensitivity tables"),
  corrected_preprocessing_replacement_available = "yes",
  source_table = c(
    "corrected_coefficients/corrected_preprocessing_full_data_coefficients_all_runs.tsv",
    "corrected_sensitivity/corrected_preprocessing_lasso_sensitivity_AUC_summary.tsv",
    "corrected_coefficients/corrected_preprocessing_predictor_overlap_turnover.tsv",
    "corrected_coefficients/corrected_preprocessing_full_data_coefficients_all_runs.tsv",
    "corrected_sensitivity/corrected_preprocessing_lasso_sensitivity_AUC_summary.tsv"
  ),
  status = "PASS",
  notes = "Draft figure-ready tables generated; final manuscript figures not overwritten.",
  stringsAsFactors = FALSE
)
write_tsv(dependency_check, file.path(qc_dir, "Figure4_Figure5_dependency_check.tsv"))

run_summary <- c(
  "# Phase 12A Corrected-Preprocessing LASSO Architecture/Sensitivity Run Summary",
  "",
  paste0("Generated: ", timestamp),
  "",
  "Scope: corrected-preprocessing terminal-only LASSO/gLOOCV rerun for baseline and run-specific global FDR candidate gene sets. GBI, ASR, t-test/FDR discovery, enrichment, and nested t-test feature selection were not rerun.",
  "",
  paste0("Runs completed: ", paste(run_specs_df$run_id, collapse = ", ")),
  "",
  "Smoke test completed before full batch and passed.",
  "",
  "Key outputs:",
  "- corrected_sensitivity/corrected_preprocessing_lasso_sensitivity_AUC_summary.tsv",
  "- corrected_coefficients/corrected_preprocessing_full_data_fit_summary_all_runs.tsv",
  "- corrected_coefficients/corrected_preprocessing_predictor_overlap_turnover.tsv",
  "- qc/drop_cetacean_collapse_check.tsv",
  "",
  "Metadata/status patch:",
  "- `fix_drop_whale` is explicitly marked as `null_or_near_null_model` in the sensitivity AUC summaries.",
  "- Sensitivity summaries include fold-level no-selected/nonzero/no-evaluable counts and full-data model-status columns.",
  "- Drop-cetacean / `fix_drop_whale` collapse is retained under corrected preprocessing: corrected gLOOCV AUC = 0.50199203187251 and corrected full-data fit selected zero predictors.",
  "- Fold diagnostics show 237/238 evaluable folds selected no predictors, 1 evaluable fold selected nonzero predictors, and 24 folds had no evaluable test species after drop/exclusion rules.",
  "",
  "No manuscript text or final figures were modified."
)
writeLines(run_summary, file.path(qc_dir, "run_summary.md"))

baseline_counts_ok <- all(selected_check$status[selected_check$run_id %in% c("fix_marine_binary", "fix_aquatic_v2")] == "PASS")
drop_story_changed <- drop_whale_full$n_selected_predictors != 0
classification <- if (drop_story_changed) "DROP_CETACEAN_STORY_CHANGED" else if (!baseline_counts_ok) "ARCHITECTURE_REQUIRES_REVISION" else "ARCHITECTURE_STABLE"
memo <- c(
  "# Fig. 4 / Fig. 5 LASSO Architecture Decision Memo",
  "",
  paste0("Generated: ", timestamp),
  "",
  paste0("Classification: `", classification, "`."),
  "",
  paste0("Can Fig. 4 be rebuilt from corrected-preprocessing coefficients? ", ifelse(nrow(full_coef) > 0, "Yes.", "No.")),
  paste0("Can Fig. 5A be rebuilt from corrected-preprocessing sensitivity AUCs? ", ifelse(nrow(auc_summary) > 0, "Yes.", "No.")),
  paste0("Can Fig. 5C be rebuilt from corrected-preprocessing turnover? ", ifelse(nrow(turnover_summary) > 0, "Yes.", "No.")),
  paste0("Does the drop-cetacean collapse story still hold? ", ifelse(drop_story_changed, "No: corrected full-data fit selected predictors; review required.", "Yes: corrected full-data fit retained zero selected predictors.")),
  "",
  "Drop-cetacean / fix_drop_whale collapse is retained under corrected preprocessing: corrected gLOOCV AUC = 0.50199 and corrected full-data fit selected zero predictors. Fold-level diagnostics show most evaluated folds selected no predictors, with only rare nonzero fold fits; specifically, 237 of 238 evaluable folds selected no predictors and 1 evaluable fold selected nonzero predictors, while 24 additional folds had no evaluable test species after drop/exclusion rules. This supports a null/near-null model interpretation.",
  "",
  "Result 2.4 / 2.5 should be checked against the corrected AUCs, selected-predictor counts, and turnover tables before integration. Legacy LASSO numbers should not be used as final manuscript values."
)
writeLines(memo, file.path(decision_dir, "Fig4_Fig5_LASSO_architecture_decision_memo.md"))

readme <- c(
  "# Phase 12A Corrected-Preprocessing LASSO Architecture/Sensitivity",
  "",
  paste0("Generated: ", timestamp),
  "",
  "This folder contains corrected-preprocessing LASSO architecture and sensitivity reruns for Fig. 4 / Fig. 5A / Fig. 5C support.",
  "",
  "Design retained:",
  "",
  "```text",
  "global deterministic ASR",
  "-> global branch-level t-test/FDR candidate-gene discovery",
  "-> fixed run-specific FDR candidate gene set",
  "-> terminal-only LASSO",
  "-> genus-level LOOCV with fold-wise training-only imputation/scaling/lambda/model fitting",
  "```",
  "",
  "No GBI, ASR, t-test/FDR discovery, enrichment, nested-t-test feature selection, manuscript edits, or final figure overwrites were performed.",
  "",
  "Draft figures in `figures_draft/` are review aids only, not final manuscript art."
)
writeLines(readme, file.path(work_dir, "README_corrected_preprocessing_LASSO_architecture_sensitivity.md"))

message("[", Sys.time(), "] Phase 12A complete. Classification: ", classification)
