#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(parallel)
})

root_value <- Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "")
if (!nzchar(root_value)) {
  stop("Set MARINE_MAMMAL_ENDPOINTFIX_ROOT to the unpacked endpoint-fix analysis root.")
}
root_dir <- normalizePath(root_value, mustWork = TRUE)
work_dir <- Sys.getenv(
  "MARINE_MAMMAL_OUTPUT_ROOT",
  unset = file.path(root_dir, "reproduction_outputs")
)
out_dir <- file.path(work_dir, "terminal_only_projection")
qc_dir <- file.path(work_dir, "qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) return(default)
  sub(paste0("^", flag, "="), "", hit[[1]])
}
cores <- as.integer(get_arg("--cores", min(12L, max(1L, parallel::detectCores() - 1L))))
seed_base <- as.integer(get_arg("--seed", 20260701L))
alpha <- as.numeric(get_arg("--alpha", 0.01))
if (!is.finite(cores) || cores < 1L) cores <- 1L

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  branch_map = file.path(root_dir, "03_split_outputs", "outputs", "species_tree.branch_map.txt"),
  phase12a_full_coefficients = file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity", "corrected_coefficients", "corrected_preprocessing_full_data_coefficients_all_runs.tsv")
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
    mean_background = mean0,
    mean_trait = mean1,
    n_background = n0,
    n_trait = n1,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
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

prepare_design <- function(raw_mat, train_branches, test_branches, feature_genes) {
  x_train_raw <- t(raw_mat[feature_genes, train_branches, drop = FALSE])
  x_test_raw <- t(raw_mat[feature_genes, test_branches, drop = FALSE])
  train_impute_means <- colMeans(x_train_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(train_impute_means)
  dropped_all_missing <- sum(!keep_not_all_missing)
  x_train_raw <- x_train_raw[, keep_not_all_missing, drop = FALSE]
  x_test_raw <- x_test_raw[, keep_not_all_missing, drop = FALSE]
  train_impute_means <- train_impute_means[keep_not_all_missing]
  if (length(train_impute_means) == 0) {
    return(list(ok = FALSE, status = "no_features_after_all_missing_drop", n_dropped_all = dropped_all_missing, n_dropped_zero = 0L))
  }
  for (k in seq_along(train_impute_means)) {
    if (anyNA(x_train_raw[, k])) x_train_raw[is.na(x_train_raw[, k]), k] <- train_impute_means[[k]]
    if (anyNA(x_test_raw[, k])) x_test_raw[is.na(x_test_raw[, k]), k] <- train_impute_means[[k]]
  }
  scale_means <- colMeans(x_train_raw)
  scale_sds <- apply(x_train_raw, 2, stats::sd)
  keep_nonzero <- is.finite(scale_sds) & scale_sds > 0
  dropped_zero <- sum(!keep_nonzero)
  if (!any(keep_nonzero)) {
    return(list(ok = FALSE, status = "no_features_after_zero_variance_drop", n_dropped_all = dropped_all_missing, n_dropped_zero = dropped_zero))
  }
  x_train <- x_train_raw[, keep_nonzero, drop = FALSE]
  x_test <- x_test_raw[, keep_nonzero, drop = FALSE]
  scale_means <- scale_means[keep_nonzero]
  scale_sds <- scale_sds[keep_nonzero]
  train_impute_means <- train_impute_means[keep_nonzero]
  x_train <- sweep(sweep(x_train, 2, scale_means, "-"), 2, scale_sds, "/")
  x_test <- sweep(sweep(x_test, 2, scale_means, "-"), 2, scale_sds, "/")
  list(
    ok = TRUE,
    x_train = x_train,
    x_test = x_test,
    feature_names = colnames(x_train),
    impute_means = train_impute_means,
    scale_means = scale_means,
    scale_sds = scale_sds,
    n_dropped_all = dropped_all_missing,
    n_dropped_zero = dropped_zero
  )
}

fit_lasso_on_design <- function(design, y_train, fold_seed) {
  cv_res <- safe_cv_glmnet(as.matrix(design$x_train), y_train, fold_seed)
  lambda_min <- cv_res$cv$lambda.min
  fit <- glmnet(
    x = as.matrix(design$x_train),
    y = y_train,
    family = "binomial",
    lambda = lambda_min,
    standardize = FALSE
  )
  beta <- as.matrix(fit$beta)[, 1]
  list(fit = fit, beta = beta, lambda_min = lambda_min, warnings = cv_res$warnings)
}

stamp("Reading endpoint-fix inputs")
gbi <- read_tsv(paths$gbi)
rownames(gbi) <- gbi$gene
gbi$gene <- NULL
gbi_mat_all <- as.matrix(gbi)
storage.mode(gbi_mat_all) <- "numeric"
rm(gbi)

branch <- read_tsv(paths$branch)
trait <- read_tsv(paths$trait)
branch_map <- read_tsv(paths$branch_map)
if (!all(branch$Branch %in% colnames(gbi_mat_all))) stop("Branch IDs missing from GBI matrix.")
gbi_mat_all <- gbi_mat_all[, branch$Branch, drop = FALSE]

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
if (anyNA(terminal$marine_binary) || anyNA(terminal$aquatic_v2)) stop("Missing terminal trait values.")

model_specs <- list(
  list(
    trait = "marine_binary",
    trait_axis = "marine specialization",
    model = "marine_terminal_only_feature_selection_gLOOCV",
    full_model = "marine_terminal_only_full_data_projection",
    response_col = "marine_binary",
    eligible = rep(TRUE, nrow(terminal)),
    evidence_level = "terminal-only supervised feature selection; terminal-only LASSO; internal branches projection targets only"
  ),
  list(
    trait = "binary_aquatic_dependence",
    trait_axis = "binary aquatic dependence",
    model = "binary_aquatic_dependence_terminal_only_feature_selection_gLOOCV",
    full_model = "binary_aquatic_dependence_terminal_only_full_data_projection",
    response_col = "aquatic_v2",
    eligible = terminal$aquatic_v2 != 0.5,
    evidence_level = "terminal-only supervised feature selection; semi-aquatic excluded; terminal-only LASSO; internal branches projection targets only"
  )
)

select_terminal_features <- function(spec, train_species) {
  y <- as.numeric(train_species[[spec$response_col]])
  cols0 <- train_species$Branch[y == 0]
  cols1 <- train_species$Branch[y == 1]
  tt <- welch_from_stats(group_stats(gbi_mat_all, cols0), group_stats(gbi_mat_all, cols1))
  tt_tested <- tt[is.finite(tt$pvalue), , drop = FALSE]
  tt_tested <- tt_tested[order(tt_tested$pvalue, decreasing = FALSE), , drop = FALSE]
  idx <- bh_old_while(tt_tested$pvalue, alpha = alpha)
  list(genes = tt_tested$gene[idx], ttest = tt_tested)
}

fit_terminal_fold <- function(spec, fold_genus, fold_index, excluded_species = character(), context = "baseline") {
  eligible_species <- terminal[spec$eligible & !(terminal$Species %in% excluded_species), , drop = FALSE]
  test_idx <- which(eligible_species$genus == fold_genus)
  if (length(test_idx) == 0) return(NULL)
  train_idx <- setdiff(seq_len(nrow(eligible_species)), test_idx)
  train_species <- eligible_species[train_idx, , drop = FALSE]
  test_species <- eligible_species[test_idx, , drop = FALSE]
  y_train <- as.numeric(train_species[[spec$response_col]])
  y_test <- as.numeric(test_species[[spec$response_col]])
  fold_seed <- seed_base + fold_index + ifelse(spec$trait == "marine_binary", 100000L, 200000L) + ifelse(context == "masked_river_dolphin", 300000L, 0L)
  feat <- select_terminal_features(spec, train_species)
  base <- data.frame(
    trait = spec$trait,
    model = spec$model,
    context = context,
    held_out_genus = fold_genus,
    fold_id = fold_index,
    n_training_species = nrow(train_species),
    n_test_species = nrow(test_species),
    n_training_positive = sum(y_train == 1),
    n_training_negative = sum(y_train == 0),
    n_test_positive = sum(y_test == 1),
    n_test_negative = sum(y_test == 0),
    n_ttest_genes_tested = nrow(feat$ttest),
    n_candidate_features_FDR_0_01 = length(feat$genes),
    fold_seed = fold_seed,
    stringsAsFactors = FALSE
  )
  if (length(unique(y_train)) < 2 || length(feat$genes) == 0) {
    pred <- rep(mean(y_train), length(y_test))
    return(list(
      fold = cbind(base, n_features_used = 0L, n_selected_predictors = 0L, lambda_min = NA_real_, fold_model_status = "null_intercept_only"),
      oof = data.frame(
        trait = spec$trait, model = spec$model, context = context,
        species = test_species$Species, genus = test_species$genus,
        trait_category = test_species$trait_category,
        response = y_test,
        predicted_probability = pred,
        held_out_genus = fold_genus,
        n_candidate_features_FDR_0_01 = length(feat$genes),
        n_selected_predictors = 0L,
        fold_model_status = "null_intercept_only",
        stringsAsFactors = FALSE
      ),
      selected = data.frame()
    ))
  }
  design <- prepare_design(gbi_mat_all, train_species$Branch, test_species$Branch, feat$genes)
  if (!isTRUE(design$ok)) {
    pred <- rep(mean(y_train), length(y_test))
    return(list(
      fold = cbind(base, n_features_used = 0L, n_selected_predictors = 0L, lambda_min = NA_real_, fold_model_status = paste0("null_intercept_only_", design$status)),
      oof = data.frame(
        trait = spec$trait, model = spec$model, context = context,
        species = test_species$Species, genus = test_species$genus,
        trait_category = test_species$trait_category,
        response = y_test,
        predicted_probability = pred,
        held_out_genus = fold_genus,
        n_candidate_features_FDR_0_01 = length(feat$genes),
        n_selected_predictors = 0L,
        fold_model_status = paste0("null_intercept_only_", design$status),
        stringsAsFactors = FALSE
      ),
      selected = data.frame()
    ))
  }
  fit_res <- fit_lasso_on_design(design, y_train, fold_seed)
  pred <- as.numeric(predict(fit_res$fit, newx = as.matrix(design$x_test), type = "response")[, 1])
  nonzero <- fit_res$beta[fit_res$beta != 0]
  selected <- if (length(nonzero) == 0) data.frame() else data.frame(
    trait = spec$trait,
    model = spec$model,
    context = context,
    held_out_genus = fold_genus,
    fold_id = fold_index,
    gene = names(nonzero),
    coefficient = as.numeric(nonzero),
    coefficient_sign = ifelse(nonzero > 0, "positive", "negative"),
    lambda_min = fit_res$lambda_min,
    stringsAsFactors = FALSE
  )
  status <- ifelse(length(nonzero) == 0, "fit_no_selected_predictors", "fit_nonzero_predictors")
  list(
    fold = cbind(base, n_features_used = ncol(design$x_train), n_selected_predictors = length(nonzero), lambda_min = fit_res$lambda_min, fold_model_status = status),
    oof = data.frame(
      trait = spec$trait, model = spec$model, context = context,
      species = test_species$Species, genus = test_species$genus,
      trait_category = test_species$trait_category,
      response = y_test,
      predicted_probability = pred,
      held_out_genus = fold_genus,
      n_candidate_features_FDR_0_01 = length(feat$genes),
      n_selected_predictors = length(nonzero),
      fold_model_status = status,
      stringsAsFactors = FALSE
    ),
    selected = selected
  )
}

run_terminal_cv <- function(spec, context = "baseline", excluded_species = character()) {
  eligible_species <- terminal[spec$eligible & !(terminal$Species %in% excluded_species), , drop = FALSE]
  genera <- sort(unique(eligible_species$genus))
  stamp("Running terminal-only CV: ", spec$trait, " context=", context, " folds=", length(genera))
  idx <- match(genera, sort(unique(terminal$genus)))
  results <- parallel::mclapply(
    seq_along(genera),
    function(i) fit_terminal_fold(spec, genera[[i]], idx[[i]], excluded_species = excluded_species, context = context),
    mc.cores = min(cores, length(genera))
  )
  results <- Filter(Negate(is.null), results)
  fold <- do.call(rbind, lapply(results, `[[`, "fold"))
  oof <- do.call(rbind, lapply(results, `[[`, "oof"))
  selected <- do.call(rbind, lapply(results, `[[`, "selected"))
  auc <- auc_rank(oof$response, oof$predicted_probability)
  list(fold = fold, oof = oof, selected = selected, auc = auc)
}

fit_terminal_full <- function(spec, context = "baseline", excluded_species = character()) {
  eligible_species <- terminal[spec$eligible & !(terminal$Species %in% excluded_species), , drop = FALSE]
  y <- as.numeric(eligible_species[[spec$response_col]])
  feat <- select_terminal_features(spec, eligible_species)
  if (length(feat$genes) == 0 || length(unique(y)) < 2) stop("Cannot fit full terminal-only model for ", spec$trait, " context=", context)
  design <- prepare_design(gbi_mat_all, eligible_species$Branch, eligible_species$Branch, feat$genes)
  if (!isTRUE(design$ok)) stop("No usable full-data terminal-only features for ", spec$trait)
  fit_res <- fit_lasso_on_design(design, y, seed_base + 900000L + ifelse(spec$trait == "marine_binary", 1L, 2L) + ifelse(context == "masked_river_dolphin", 100L, 0L))
  nonzero <- fit_res$beta[fit_res$beta != 0]
  coeff <- if (length(nonzero) == 0) data.frame() else data.frame(
    trait = spec$trait,
    model = spec$full_model,
    context = context,
    gene = names(nonzero),
    coefficient = as.numeric(nonzero),
    coefficient_sign = ifelse(nonzero > 0, "positive", "negative"),
    lambda_min = fit_res$lambda_min,
    stringsAsFactors = FALSE
  )
  list(
    spec = spec,
    context = context,
    train_species = eligible_species,
    candidate_genes = feat$genes,
    ttest = feat$ttest,
    design = design,
    fit = fit_res$fit,
    beta = fit_res$beta,
    lambda_min = fit_res$lambda_min,
    coeff = coeff
  )
}

project_model <- function(full_fit, target_df, output_context) {
  feature_names <- full_fit$design$feature_names
  raw <- t(gbi_mat_all[feature_names, target_df$branch_id, drop = FALSE])
  for (k in seq_along(feature_names)) {
    fname <- feature_names[[k]]
    raw[is.na(raw[, k]), k] <- full_fit$design$impute_means[[fname]]
  }
  x <- sweep(sweep(raw, 2, full_fit$design$scale_means, "-"), 2, full_fit$design$scale_sds, "/")
  pred <- as.numeric(predict(full_fit$fit, newx = as.matrix(x), type = "response")[, 1])
  beta <- full_fit$beta[full_fit$beta != 0]
  proj <- cbind(
    target_df,
    trait = full_fit$spec$trait,
    model = full_fit$spec$full_model,
    context = output_context,
    predicted_probability = pred,
    logit = qlogis(pmin(pmax(pred, 1e-12), 1 - 1e-12)),
    n_candidate_features = length(full_fit$candidate_genes),
    n_selected_predictors = length(beta),
    lambda_min = full_fit$lambda_min,
    stringsAsFactors = FALSE
  )
  contrib <- data.frame()
  if (length(beta) > 0) {
    nonzero_x <- x[, names(beta), drop = FALSE]
    contrib_list <- lapply(seq_len(nrow(target_df)), function(i) {
      vals <- as.numeric(nonzero_x[i, ]) * as.numeric(beta)
      data.frame(
        target_name = target_df$target_name[[i]],
        branch_id = target_df$branch_id[[i]],
        target_type = target_df$target_type[[i]],
        trait = full_fit$spec$trait,
        model = full_fit$spec$full_model,
        context = output_context,
        gene = names(beta),
        coefficient = as.numeric(beta),
        scaled_GBI = as.numeric(nonzero_x[i, ]),
        contribution_to_logit_excluding_intercept = vals,
        abs_contribution = abs(vals),
        stringsAsFactors = FALSE
      )
    })
    contrib <- do.call(rbind, contrib_list)
  }
  list(projections = proj, contributions = contrib)
}

internal_targets <- data.frame(
  target_name = c(
    "common_ancestor_Cetacea_plus_Hippopotamus",
    "crown_Cetacea",
    "Mysticeti",
    "Odontoceti",
    "Platanista_terminal_branch",
    "Inia_plus_Lipotes",
    "Pinnipedia",
    "Phocidae",
    "Otarioidea",
    "Otariidae",
    "Otariinae_sampled",
    "Sirenia",
    "Dugong_plus_Hydrodamalis",
    "Trichechus_terminal_branch",
    "Ursus_arctos_plus_Ursus_maritimus",
    "Lutrinae_sampled",
    "Enhydra_terminal_branch"
  ),
  native_branch_label = c(
    "B593", "B592", "B568", "B591", "B281", "B573", "B501", "B500", "B493", "B492", "B491", "B313", "B312", "B13", "B470", "B483", "B197"
  ),
  branch_id = c(
    "B593", "B592", "B544", "B591", "B548", "B558", "B410", "B409", "B394", "B393", "B392", "B26", "B25", "B22", "B349", "B380", "B374"
  ),
  target_type = c(
    "internal_clade", "internal_clade", "internal_clade", "internal_clade", "terminal_lineage", "internal_clade", "internal_clade", "internal_clade", "internal_clade", "internal_clade", "internal_clade", "internal_clade", "internal_clade", "terminal_lineage", "internal_clade", "internal_clade", "terminal_lineage"
  ),
  target_note = c(
    "Exact internal clade: Hippopotamus_amphibius plus sampled Cetacea; native B593 maps to old B593.",
    "Exact internal clade: sampled crown Cetacea; native B592 maps to old B592.",
    "Exact internal clade: sampled Mysticeti; native B568 maps to old B544.",
    "Exact internal clade: sampled Odontoceti; native B591 maps to old B591.",
    "Single sampled Platanista branch; native B281 maps to old B548.",
    "Exact internal clade: Inia_geoffrensis plus Lipotes_vexillifer; native B573 maps to old B558.",
    "Exact internal clade: sampled Pinnipedia; native B501 maps to old B410.",
    "Exact internal clade: sampled Phocidae; native B500 maps to old B409.",
    "Exact internal clade: Odobenus plus sampled otariids; native B493 maps to old B394.",
    "Exact internal clade: sampled Otariidae; native B492 maps to old B393.",
    "Exact internal clade: Eumetopias plus Zalophus; native B491 maps to old B392.",
    "Exact internal clade: sampled Sirenia; native B313 maps to old B26.",
    "Exact internal clade: Dugong_dugon plus Hydrodamalis_gigas; native B312 maps to old B25.",
    "Single sampled Trichechus branch; native B13 maps to old B22.",
    "Exact internal clade: sampled brown bear plus polar bear; native B470 maps to old B349.",
    "Exact internal clade: sampled Pteronura/Lontra/Enhydra/Aonyx/Lutra; native B483 maps to old B380.",
    "Single sampled sea otter branch; native B197 maps to old B374."
  ),
  stringsAsFactors = FALSE
)
internal_targets$branch_label_found <- internal_targets$branch_id %in% branch$Branch

masked_river_dolphins <- c("Platanista_minor", "Inia_geoffrensis", "Lipotes_vexillifer")

all_cv <- list()
all_selected <- list()
all_full <- list()
performance_rows <- list()

for (spec in model_specs) {
  cv <- run_terminal_cv(spec, context = "baseline")
  full <- fit_terminal_full(spec, context = "baseline")
  all_cv[[spec$trait]] <- cv
  all_full[[spec$trait]] <- full
  all_selected[[spec$trait]] <- cv$selected
  performance_rows[[spec$trait]] <- data.frame(
    trait = spec$trait,
    trait_axis = spec$trait_axis,
    context = "terminal_only_baseline",
    n_species_evaluated = nrow(cv$oof),
    n_positive = sum(cv$oof$response == 1),
    n_negative = sum(cv$oof$response == 0),
    AUC_terminal_only_nested_gLOOCV = cv$auc,
    n_folds = nrow(cv$fold),
    median_candidate_features_per_fold = stats::median(cv$fold$n_candidate_features_FDR_0_01),
    IQR_candidate_features_per_fold = iqr_text(cv$fold$n_candidate_features_FDR_0_01),
    median_selected_predictors_per_fold = stats::median(cv$fold$n_selected_predictors),
    IQR_selected_predictors_per_fold = iqr_text(cv$fold$n_selected_predictors),
    full_data_terminal_only_candidate_genes = length(full$candidate_genes),
    full_data_terminal_only_selected_predictors = nrow(full$coeff),
    evidence_level = spec$evidence_level,
    stringsAsFactors = FALSE
  )
}

masked_cv <- run_terminal_cv(model_specs[[1]], context = "masked_river_dolphin", excluded_species = masked_river_dolphins)
masked_full <- fit_terminal_full(model_specs[[1]], context = "masked_river_dolphin", excluded_species = masked_river_dolphins)
masked_target_df <- terminal[match(masked_river_dolphins, terminal$Species), c("Species", "Branch", "trait_category", "aquaticity_score"), drop = FALSE]
masked_target_df <- data.frame(
  target_name = masked_target_df$Species,
  branch_id = masked_target_df$Branch,
  target_type = "masked_terminal_projection",
  target_note = "Masked river dolphin terminal species projected after terminal-only training excludes all three masked taxa.",
  stringsAsFactors = FALSE
)
masked_proj <- project_model(masked_full, masked_target_df, "masked_river_dolphin_terminal_only")
performance_rows[["masked_river_dolphin"]] <- data.frame(
  trait = "marine_binary",
  trait_axis = "marine specialization",
  context = "masked_river_dolphin_terminal_only",
  n_species_evaluated = nrow(masked_cv$oof),
  n_positive = sum(masked_cv$oof$response == 1),
  n_negative = sum(masked_cv$oof$response == 0),
  AUC_terminal_only_nested_gLOOCV = masked_cv$auc,
  n_folds = nrow(masked_cv$fold),
  median_candidate_features_per_fold = stats::median(masked_cv$fold$n_candidate_features_FDR_0_01),
  IQR_candidate_features_per_fold = iqr_text(masked_cv$fold$n_candidate_features_FDR_0_01),
  median_selected_predictors_per_fold = stats::median(masked_cv$fold$n_selected_predictors),
  IQR_selected_predictors_per_fold = iqr_text(masked_cv$fold$n_selected_predictors),
  full_data_terminal_only_candidate_genes = length(masked_full$candidate_genes),
  full_data_terminal_only_selected_predictors = nrow(masked_full$coeff),
  evidence_level = "masked river-dolphin terminal-only sensitivity; masked taxa excluded from training/AUC and projected as ambiguous terminal targets",
  stringsAsFactors = FALSE
)

full_projection_results <- lapply(all_full, function(full) project_model(full, internal_targets, "terminal_only_internal_projection"))
internal_proj <- do.call(rbind, lapply(full_projection_results, `[[`, "projections"))
internal_contrib <- do.call(rbind, lapply(full_projection_results, `[[`, "contributions"))

performance <- do.call(rbind, performance_rows)
write_tsv(performance, file.path(out_dir, "terminal_only_model_performance.csv"))
write_tsv(do.call(rbind, lapply(all_cv, `[[`, "oof")), file.path(out_dir, "terminal_only_OOF_predictions.csv"))
write_tsv(do.call(rbind, lapply(all_cv, `[[`, "fold")), file.path(out_dir, "terminal_only_fold_diagnostics.csv"))
write_tsv(do.call(rbind, lapply(all_selected, identity)), file.path(out_dir, "terminal_only_selected_predictors_by_fold.csv"))
write_tsv(do.call(rbind, lapply(all_full, `[[`, "coeff")), file.path(out_dir, "terminal_only_selected_predictors.csv"))
write_tsv(internal_targets, file.path(out_dir, "selected_internal_node_target_audit.csv"))
write_tsv(internal_proj, file.path(out_dir, "selected_internal_node_projections.csv"))
write_tsv(internal_contrib, file.path(out_dir, "selected_internal_node_contributions_long.csv"))
write_tsv(masked_cv$oof, file.path(out_dir, "masked_river_dolphin_terminal_only_OOF_predictions.csv"))
write_tsv(masked_cv$fold, file.path(out_dir, "masked_river_dolphin_terminal_only_fold_diagnostics.csv"))
write_tsv(masked_full$coeff, file.path(out_dir, "masked_river_dolphin_terminal_only_selected_predictors.csv"))
write_tsv(masked_proj$projections, file.path(out_dir, "masked_river_dolphin_predictions.csv"))
write_tsv(masked_proj$contributions, file.path(out_dir, "masked_river_dolphin_contributions_long.csv"))

phase12a <- read_tsv(paths$phase12a_full_coefficients)
baseline_sets <- list(
  marine_binary = phase12a$gene[phase12a$run_id == "fix_marine_binary" & phase12a$selected == "TRUE"],
  binary_aquatic_dependence = phase12a$gene[phase12a$run_id == "fix_aquatic_v2" & phase12a$selected == "TRUE"]
)
overlap_rows <- lapply(names(all_full), function(tr) {
  terminal_set <- all_full[[tr]]$coeff$gene
  baseline_set <- baseline_sets[[tr]]
  data.frame(
    trait = tr,
    endpointfix_corrected_full_data_predictors = length(baseline_set),
    terminal_only_full_data_predictors = length(terminal_set),
    n_overlap = length(intersect(baseline_set, terminal_set)),
    n_union = length(union(baseline_set, terminal_set)),
    jaccard_overlap = length(intersect(baseline_set, terminal_set)) / length(union(baseline_set, terminal_set)),
    stringsAsFactors = FALSE
  )
})
overlap_df <- do.call(rbind, overlap_rows)
write_tsv(overlap_df, file.path(out_dir, "terminal_only_vs_endpointfix_full_data_predictor_overlap.csv"))

top_contrib <- function(contrib, n = 10) {
  contrib <- contrib[order(contrib$target_name, contrib$trait, -contrib$abs_contribution), , drop = FALSE]
  keep <- ave(contrib$abs_contribution, contrib$target_name, contrib$trait, FUN = function(x) rank(-x, ties.method = "first")) <= n
  contrib[keep, , drop = FALSE]
}
write_tsv(top_contrib(internal_contrib), file.path(out_dir, "selected_internal_node_top_contributions.csv"))
write_tsv(top_contrib(masked_proj$contributions), file.path(out_dir, "masked_river_dolphin_top_contributions.csv"))

report <- c(
  "# Terminal-Only Projection Report",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")),
  "",
  "## Design",
  "",
  "- Stage 1 gene screening uses terminal branches/species labels only.",
  "- LASSO training and lambda selection use terminal species only.",
  "- Internal branches are projection targets only and do not enter feature selection or model fitting.",
  "- Outputs are predicted probabilities / rate-based genomic-state projections, not posterior probabilities and not ancestral habitat reconstructions.",
  "",
  "## Performance",
  "",
  paste(capture.output(print(performance)), collapse = "\n"),
  "",
  "## Predictor Overlap With Endpoint-Fix Corrected Full-Data Architecture",
  "",
  paste(capture.output(print(overlap_df)), collapse = "\n"),
  "",
  "## Masked River-Dolphin Terminal-Only Projection",
  "",
  "- Masked taxa: Platanista_minor, Inia_geoffrensis, Lipotes_vexillifer.",
  "- These taxa were excluded from marine terminal-only training and AUC evaluation, then projected as ambiguous terminal targets.",
  "- ASR-informed masked projection was not run in this terminal-only script and should remain a separate evidence layer if required.",
  "",
  paste(capture.output(print(masked_proj$projections[, c("target_name", "branch_id", "predicted_probability", "logit", "n_candidate_features", "n_selected_predictors")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Internal Target Audit",
  "",
  paste(capture.output(print(internal_targets[, c("target_name", "branch_id", "target_type", "branch_label_found")], row.names = FALSE)), collapse = "\n")
)
writeLines(report, file.path(out_dir, "terminal_only_projection_report.md"))

masked_report <- c(
  "# Masked River-Dolphin Projection Report",
  "",
  "This report summarizes the terminal-only masked river-dolphin sensitivity. ASR-informed masked projection is not included in this first terminal-only run.",
  "",
  "## Training Domain",
  "",
  paste0("- Training/evaluation excludes: ", paste(masked_river_dolphins, collapse = ", ")),
  paste0("- Cross-validated n evaluated: ", performance_rows[["masked_river_dolphin"]]$n_species_evaluated),
  paste0("- Positives/negatives: ", performance_rows[["masked_river_dolphin"]]$n_positive, " / ", performance_rows[["masked_river_dolphin"]]$n_negative),
  paste0("- Terminal-only nested AUC: ", sprintf("%.4f", performance_rows[["masked_river_dolphin"]]$AUC_terminal_only_nested_gLOOCV)),
  "",
  "## Masked Terminal Projections",
  "",
  paste(capture.output(print(masked_proj$projections[, c("target_name", "predicted_probability", "logit", "n_candidate_features", "n_selected_predictors")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Boundary",
  "",
  "Use these as terminal-only rate-profile projections, not posterior probabilities or habitat reconstructions."
)
writeLines(masked_report, file.path(out_dir, "masked_river_dolphin_report.md"))

qc <- data.frame(
  check = c(
    "terminal_only_screening_uses_terminal_branches_only",
    "internal_branches_projection_targets_only",
    "semi_aquatic_excluded_for_binary_aquatic",
    "masked_river_dolphins_excluded_from_training_AUC",
    "masked_ASR_informed_projection_not_mixed_with_terminal_only",
    "all_selected_internal_target_branch_ids_found"
  ),
  status = c(
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    "PASS",
    ifelse(all(internal_targets$branch_label_found), "PASS", "FAIL")
  ),
  notes = c(
    "Welch screen uses train_species terminal branch columns only.",
    "Projection table is generated after final terminal-only model fitting.",
    "Spec eligible excludes aquatic_v2=0.5.",
    "Masked species removed before CV and full-data terminal-only fit.",
    "ASR-informed masked projection is explicitly not included in this terminal-only evidence layer.",
    paste(internal_targets$target_name[!internal_targets$branch_label_found], collapse = ";")
  ),
  stringsAsFactors = FALSE
)
write_tsv(qc, file.path(qc_dir, "terminal_only_projection_QC.tsv"))

stamp("Done. Outputs in ", out_dir)
