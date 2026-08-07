#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})

root_value <- Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "")
if (!nzchar(root_value)) {
  stop("Set MARINE_MAMMAL_ENDPOINTFIX_ROOT to the unpacked endpoint-fix analysis root.")
}
root_dir <- normalizePath(root_value, mustWork = TRUE)
output_root <- Sys.getenv(
  "MARINE_MAMMAL_NC_OUTPUT_ROOT",
  unset = file.path(root_dir, "analysis_nc_arc")
)
out_dir <- file.path(output_root, "figureS3_ancestor_fingerprints")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(file, ...) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_csv <- function(x, file) utils::write.csv(x, file = file, quote = TRUE, row.names = FALSE, na = "")
write_tsv <- function(x, file) utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
inv_logit <- function(x) 1 / (1 + exp(-x))
fmt3 <- function(x) ifelse(is.na(x), "NA", sprintf("%.3f", x))

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  marine_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_marine_binary", "marine_binary.mammal.FDR0.01.n1559.t_test.txt"),
  aquatic_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_aquatic_v2", "aquatic_v2.mammal.FDR0.01.n1227.t_test.txt"),
  archived_coef = file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity", "corrected_coefficients", "corrected_preprocessing_full_data_coefficients_all_runs.tsv"),
  table_s5 = file.path(root_dir, "supplementary_tables_endpointfix", "TableS5", "Supplementary_Table_S5_Figure4C_predictor_annotation.tsv"),
  fig6a_profiles = Sys.getenv(
    "MARINE_MAMMAL_FIG6A_PROFILE_SOURCE",
    unset = file.path(root_dir, "analysis_nc_arc", "figure4_rebalanced_main", "SourceData_Fig4A_paired_projection_profiles.csv")
  ),
  terminal_only_internal = file.path(root_dir, "analysis_nc_arc", "terminal_only_projection", "selected_internal_node_projections.csv"),
  terminal_only_audit = file.path(root_dir, "analysis_nc_arc", "terminal_only_projection", "selected_internal_node_target_audit.csv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

message("[", Sys.time(), "] Reading GBI matrix and frozen model inputs")
gbi <- read_tsv(paths$gbi)
stopifnot("gene" %in% names(gbi))
rownames(gbi) <- gbi$gene
gbi$gene <- NULL
gbi_mat <- as.matrix(gbi)
storage.mode(gbi_mat) <- "numeric"

branch <- read_tsv(paths$branch)
trait <- read_tsv(paths$trait)
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
stopifnot(all(terminal$Branch %in% colnames(gbi_mat)))

archived_coef <- read_tsv(paths$archived_coef)
table_s5 <- read_tsv(paths$table_s5)
fig6a <- utils::read.csv(paths$fig6a_profiles, stringsAsFactors = FALSE, check.names = FALSE)
terminal_only_internal <- read_tsv(paths$terminal_only_internal)
target_audit <- read_tsv(paths$terminal_only_audit)

make_spec <- function(model) {
  if (model == "marine") {
    genes <- read_tsv(paths$marine_features)$gene
    list(
      model = "marine",
      trait = "marine_binary",
      run_id = "fix_marine_binary",
      genes = genes,
      y = as.numeric(terminal$marine_binary),
      eligible = is.finite(terminal$marine_binary)
    )
  } else {
    genes <- read_tsv(paths$aquatic_features)$gene
    y <- as.numeric(terminal$aquatic_v2)
    list(
      model = "binary_aquatic_dependence",
      trait = "binary_aquatic_dependence",
      run_id = "fix_aquatic_v2",
      genes = genes,
      y = y,
      eligible = is.finite(y) & y != 0.5
    )
  }
}

fit_full_model <- function(spec) {
  eligible_species <- terminal[spec$eligible, , drop = FALSE]
  y <- as.numeric(spec$y[spec$eligible])
  x_raw <- t(gbi_mat[spec$genes, eligible_species$Branch, drop = FALSE])
  impute_means <- colMeans(x_raw, na.rm = TRUE)
  keep_not_all_missing <- is.finite(impute_means)
  x_raw <- x_raw[, keep_not_all_missing, drop = FALSE]
  impute_means <- impute_means[keep_not_all_missing]
  for (i in seq_along(impute_means)) {
    if (anyNA(x_raw[, i])) x_raw[is.na(x_raw[, i]), i] <- impute_means[[i]]
  }
  scale_means <- colMeans(x_raw)
  scale_sds <- apply(x_raw, 2, stats::sd)
  keep_nonzero <- is.finite(scale_sds) & scale_sds > 0
  x_scaled <- x_raw[, keep_nonzero, drop = FALSE]
  scale_means <- scale_means[keep_nonzero]
  scale_sds <- scale_sds[keep_nonzero]
  impute_means <- impute_means[keep_nonzero]
  x_scaled <- sweep(sweep(x_scaled, 2, scale_means, "-"), 2, scale_sds, "/")

  set.seed(20260523L + 900000L + length(spec$genes) + ifelse(spec$model == "marine", 1L, 8L))
  cv <- cv.glmnet(
    x = as.matrix(x_scaled),
    y = y,
    family = "binomial",
    foldid = seq_along(y),
    standardize = FALSE,
    type.measure = "deviance"
  )
  fit <- glmnet(
    x = as.matrix(x_scaled),
    y = y,
    family = "binomial",
    lambda = cv$lambda.min,
    standardize = FALSE
  )
  beta <- as.matrix(fit$beta)[, 1]
  coef_all <- as.matrix(stats::coef(fit))[, 1]
  selected <- names(beta)[beta != 0]
  list(
    spec = spec,
    eligible_species = eligible_species,
    impute_means = impute_means,
    scale_means = scale_means,
    scale_sds = scale_sds,
    lambda_min = cv$lambda.min,
    intercept = unname(coef_all["(Intercept)"]),
    beta = beta,
    selected = selected,
    dropped_all_missing = sum(!keep_not_all_missing),
    dropped_zero_variance = sum(!keep_nonzero)
  )
}

add_annotation <- function(df) {
  predictor_cols <- c(
    "gene", "lasso_group", "display_module", "pro_recommended_module",
    "recommended_submodule_or_function", "annotation_confidence",
    "Figure4C_display_status"
  )
  predictor_cols <- intersect(predictor_cols, names(table_s5))
  ann <- table_s5[, predictor_cols, drop = FALSE]
  names(ann)[names(ann) == "lasso_group"] <- "predictor_class"
  if (!"display_module" %in% names(ann)) ann$display_module <- NA_character_
  if (!"predictor_class" %in% names(ann)) ann$predictor_class <- NA_character_
  merge(df, ann, by = "gene", all.x = TRUE, sort = FALSE)
}

project_branch <- function(model_fit, target_row) {
  genes <- model_fit$selected
  branch_id <- target_row$branch_id
  raw <- as.numeric(gbi_mat[genes, branch_id])
  names(raw) <- genes
  imputed_flag <- is.na(raw)
  filled <- raw
  filled[imputed_flag] <- model_fit$impute_means[genes][imputed_flag]
  scaled <- (filled - model_fit$scale_means[genes]) / model_fit$scale_sds[genes]
  coef <- model_fit$beta[genes]
  contribution <- scaled * coef
  logit <- model_fit$intercept + sum(contribution)
  prob <- inv_logit(logit)
  data.frame(
    target_name = target_row$target_name,
    display_label = target_row$display_label,
    branch_id = branch_id,
    native_branch_label = target_row$native_branch_label,
    target_type = target_row$target_type,
    inclusion_set = target_row$inclusion_set,
    branch_group = target_row$branch_group,
    model = model_fit$spec$model,
    trait = model_fit$spec$trait,
    gene = genes,
    coefficient = as.numeric(coef),
    scaled_GBI = as.numeric(scaled),
    contribution = as.numeric(contribution),
    abs_contribution = abs(as.numeric(contribution)),
    contribution_sign = ifelse(contribution > 0, "positive", ifelse(contribution < 0, "negative", "zero")),
    GBI_raw_if_available = raw,
    imputed_flag = imputed_flag,
    fitted_logit = logit,
    fitted_probability = prob,
    intercept_included_in_profile_score = model_fit$intercept,
    stringsAsFactors = FALSE
  )
}

message("[", Sys.time(), "] Rebuilding Figure 6B/C corrected full-data baseline models")
models <- list(
  marine = fit_full_model(make_spec("marine")),
  binary_aquatic_dependence = fit_full_model(make_spec("binary_aquatic_dependence"))
)

coef_repro <- do.call(rbind, lapply(models, function(m) {
  arch <- archived_coef[archived_coef$run_id == m$spec$run_id, , drop = FALSE]
  common <- intersect(arch$gene, names(m$beta))
  arch_beta <- arch$coefficient[match(common, arch$gene)]
  new_beta <- m$beta[common]
  data.frame(
    model = m$spec$model,
    run_id = m$spec$run_id,
    archived_n_rows = nrow(arch),
    rebuilt_n_features_after_filter = length(m$beta),
    archived_n_selected = sum(arch$selected),
    rebuilt_n_selected = length(m$selected),
    selected_gene_set_identical = identical(sort(arch$gene[arch$selected]), sort(m$selected)),
    max_abs_coefficient_difference = max(abs(arch_beta - new_beta), na.rm = TRUE),
    lambda_archived = unique(arch$lambda_used)[1],
    lambda_rebuilt = m$lambda_min,
    intercept_rebuilt = m$intercept,
    status = ifelse(
      identical(sort(arch$gene[arch$selected]), sort(m$selected)) &&
        max(abs(arch_beta - new_beta), na.rm = TRUE) < 1e-10,
      "PASS", "STOP"
    ),
    stringsAsFactors = FALSE
  )
}))
write_tsv(coef_repro, file.path(out_dir, "ancestor_fullData_model_rebuild_coefficient_reproducibility.tsv"))
if (any(coef_repro$status != "PASS")) stop("Full-data model rebuild failed coefficient reproduction.")

node_defs <- data.frame(
  target_name = c(
    "common_ancestor_Cetacea_plus_Hippopotamus", "crown_Cetacea", "Mysticeti", "Odontoceti",
    "Pinnipedia", "Phocidae", "Otarioidea", "Otariidae", "Otariinae_sampled",
    "Sirenia", "Dugong_plus_Hydrodamalis",
    "Ursus_arctos_plus_Ursus_maritimus", "Lutrinae_sampled"
  ),
  display_label = c(
    "Cetacea + hippo ancestor", "Crown Cetacea", "Mysticeti", "Odontoceti",
    "Crown Pinnipedia", "Phocidae", "Otarioidea", "Otariidae", "Otariinae",
    "Sirenia", "Dugong + Hydrodamalis",
    "Brown bear + polar bear", "Lutrinae / Enhydra-related"
  ),
  branch_group = c(
    "Cetacean axis", "Cetacean axis", "Cetacean axis", "Cetacean axis",
    "Pinniped axis", "Pinniped axis", "Pinniped axis", "Pinniped axis", "Pinniped axis",
    "Sirenian / edge context", "Sirenian / edge context",
    "Sirenian / edge context", "Sirenian / edge context"
  ),
  inclusion_set = c(
    "required", "required", "optional", "required",
    "required", "required", "required", "required", "optional",
    "required", "required",
    "optional", "optional"
  ),
  plot_order = seq_len(13),
  stringsAsFactors = FALSE
)
targets <- merge(node_defs, target_audit, by = "target_name", all.x = TRUE, sort = FALSE)
targets <- targets[order(targets$plot_order), , drop = FALSE]
targets$branch_label_found <- targets$branch_label_found %in% TRUE | targets$branch_label_found == "TRUE"
targets$branch_id_found_in_gbi <- targets$branch_id %in% colnames(gbi_mat)
if (any(!targets$branch_id_found_in_gbi)) {
  stop("Missing branch IDs in GBI matrix: ", paste(targets$target_name[!targets$branch_id_found_in_gbi], collapse = ", "))
}
write_tsv(targets, file.path(out_dir, "ancestor_target_branch_map.tsv"))

message("[", Sys.time(), "] Projecting internal branches under 71/101 full-data architecture")
fp <- do.call(rbind, c(
  lapply(seq_len(nrow(targets)), function(i) project_branch(models$marine, targets[i, ])),
  lapply(seq_len(nrow(targets)), function(i) project_branch(models$binary_aquatic_dependence, targets[i, ]))
))
fp <- add_annotation(fp)
fp$module <- ifelse(!is.na(fp$display_module) & fp$display_module != "", fp$display_module, "Table S5-only / unassigned")
fp$predictor_class <- ifelse(is.na(fp$predictor_class) | fp$predictor_class == "", "not_classified", fp$predictor_class)

gene_order_list <- lapply(models, function(m) {
  genes <- m$selected
  coefs <- m$beta[genes]
  ord <- order(coefs, decreasing = FALSE)
  tmp <- data.frame(
    model = m$spec$model,
    gene = genes[ord],
    gene_order = seq_along(genes),
    coefficient = as.numeric(coefs[ord]),
    coefficient_sign = ifelse(coefs[ord] > 0, "positive", "negative"),
    stringsAsFactors = FALSE
  )
  add_annotation(tmp)
})
gene_order <- do.call(rbind, gene_order_list)
gene_order$module <- ifelse(!is.na(gene_order$display_module) & gene_order$display_module != "", gene_order$display_module, "Table S5-only / unassigned")
gene_order$predictor_class <- ifelse(is.na(gene_order$predictor_class) | gene_order$predictor_class == "", "not_classified", gene_order$predictor_class)
gene_order_out <- gene_order[, c("model", "gene", "gene_order", "coefficient", "coefficient_sign", "predictor_class", "module"), drop = FALSE]
fp <- merge(fp, gene_order_out[, c("model", "gene", "gene_order")], by = c("model", "gene"), all.x = TRUE, sort = FALSE)
fp <- fp[order(match(fp$model, c("marine", "binary_aquatic_dependence")), match(fp$target_name, targets$target_name), fp$gene_order), , drop = FALSE]

fingerprint_out <- fp[, c(
  "target_name", "display_label", "branch_id", "native_branch_label", "target_type",
  "inclusion_set", "branch_group", "model", "trait", "gene", "gene_order",
  "predictor_class", "module", "coefficient", "scaled_GBI", "contribution",
  "abs_contribution", "contribution_sign", "GBI_raw_if_available", "imputed_flag",
  "fitted_logit", "fitted_probability", "intercept_included_in_profile_score"
), drop = FALSE]
write_csv(fingerprint_out, file.path(out_dir, "SourceData_FigS3_ancestor_fingerprints_long.csv"))

profile_rows <- do.call(rbind, lapply(split(fp, paste(fp$target_name, fp$model, sep = "||")), function(part) {
  first <- part[1, , drop = FALSE]
  data.frame(
    target_name = first$target_name,
    display_label = first$display_label,
    branch_id = first$branch_id,
    native_branch_label = first$native_branch_label,
    target_type = first$target_type,
    inclusion_set = first$inclusion_set,
    branch_group = first$branch_group,
    model = first$model,
    n_predictors = nrow(part),
    n_imputed_predictors = sum(part$imputed_flag),
    intercept = first$intercept_included_in_profile_score,
    contribution_sum = sum(part$contribution),
    reconstructed_logit = first$intercept_included_in_profile_score + sum(part$contribution),
    fitted_logit = first$fitted_logit,
    fitted_probability = first$fitted_probability,
    abs_logit_reconstruction_error = abs((first$intercept_included_in_profile_score + sum(part$contribution)) - first$fitted_logit),
    evidence_layer = "Figure6B/C corrected full-data architecture; internal branch descriptive projection",
    stringsAsFactors = FALSE
  )
}))
profile_wide <- reshape(
  profile_rows[, c("target_name", "display_label", "branch_id", "native_branch_label", "target_type", "inclusion_set", "branch_group", "model", "n_predictors", "n_imputed_predictors", "fitted_probability", "fitted_logit")],
  idvar = c("target_name", "display_label", "branch_id", "native_branch_label", "target_type", "inclusion_set", "branch_group"),
  timevar = "model",
  direction = "wide"
)
names(profile_wide) <- sub("fitted_probability.marine", "marine_profile_score", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("fitted_probability.binary_aquatic_dependence", "aquatic_profile_score", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("fitted_logit.marine", "marine_logit", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("fitted_logit.binary_aquatic_dependence", "aquatic_logit", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("n_predictors.marine", "marine_n_predictors", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("n_predictors.binary_aquatic_dependence", "aquatic_n_predictors", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("n_imputed_predictors.marine", "marine_n_imputed_predictors", names(profile_wide), fixed = TRUE)
names(profile_wide) <- sub("n_imputed_predictors.binary_aquatic_dependence", "aquatic_n_imputed_predictors", names(profile_wide), fixed = TRUE)
profile_wide <- merge(profile_wide, targets[, c("target_name", "plot_order")], by = "target_name", all.x = TRUE, sort = FALSE)
fig6a_internal <- fig6a[fig6a$row_type == "internal", c("source_id", "marine_probability", "aquatic_probability", "source_layer")]
names(fig6a_internal) <- c("target_name", "figure6A_source_marine_probability", "figure6A_source_aquatic_probability", "figure6A_source_layer")
profile_wide <- merge(profile_wide, fig6a_internal, by = "target_name", all.x = TRUE, sort = FALSE)
profile_wide$abs_delta_marine_vs_figure6A_source <- abs(profile_wide$marine_profile_score - profile_wide$figure6A_source_marine_probability)
profile_wide$abs_delta_aquatic_vs_figure6A_source <- abs(profile_wide$aquatic_profile_score - profile_wide$figure6A_source_aquatic_probability)
profile_wide$figure6A_source_comparison_status <- ifelse(
  is.na(profile_wide$figure6A_source_marine_probability),
  "not_in_current_Figure6A_source",
  ifelse(
    profile_wide$abs_delta_marine_vs_figure6A_source < 1e-6 &
      profile_wide$abs_delta_aquatic_vs_figure6A_source < 1e-6,
    "PASS_matches_current_Figure6A_source",
    "REVIEW_differs_current_Figure6A_source_due_to_71_101_architecture"
  )
)
profile_wide <- profile_wide[order(profile_wide$plot_order), , drop = FALSE]
write_csv(profile_wide, file.path(out_dir, "SourceData_FigS3_ancestor_profile_scores.csv"))

top_abs <- do.call(rbind, lapply(split(fp, paste(fp$target_name, fp$model, sep = "||")), function(part) {
  part <- part[order(part$abs_contribution, decreasing = TRUE), , drop = FALSE]
  part$top_abs_rank <- seq_len(nrow(part))
  part[part$top_abs_rank <= 12, , drop = FALSE]
}))
top_pos <- do.call(rbind, lapply(split(fp[fp$contribution > 0, ], paste(fp$target_name[fp$contribution > 0], fp$model[fp$contribution > 0], sep = "||")), function(part) {
  part <- part[order(part$contribution, decreasing = TRUE), , drop = FALSE]
  part$direction_rank <- seq_len(nrow(part))
  part[part$direction_rank <= 7, , drop = FALSE]
}))
top_neg <- do.call(rbind, lapply(split(fp[fp$contribution < 0, ], paste(fp$target_name[fp$contribution < 0], fp$model[fp$contribution < 0], sep = "||")), function(part) {
  part <- part[order(part$contribution, decreasing = FALSE), , drop = FALSE]
  part$direction_rank <- seq_len(nrow(part))
  part[part$direction_rank <= 7, , drop = FALSE]
}))
top_abs$top_set <- "top_absolute"
top_pos$top_set <- "top_positive"
top_neg$top_set <- "top_negative"
top_out <- rbind(
  top_abs[, c("target_name", "display_label", "branch_id", "model", "gene", "gene_order", "module", "predictor_class", "coefficient", "scaled_GBI", "contribution", "abs_contribution", "top_abs_rank", "top_set")],
  transform(top_pos[, c("target_name", "display_label", "branch_id", "model", "gene", "gene_order", "module", "predictor_class", "coefficient", "scaled_GBI", "contribution", "abs_contribution", "direction_rank", "top_set")], top_abs_rank = direction_rank)[, c("target_name", "display_label", "branch_id", "model", "gene", "gene_order", "module", "predictor_class", "coefficient", "scaled_GBI", "contribution", "abs_contribution", "top_abs_rank", "top_set")],
  transform(top_neg[, c("target_name", "display_label", "branch_id", "model", "gene", "gene_order", "module", "predictor_class", "coefficient", "scaled_GBI", "contribution", "abs_contribution", "direction_rank", "top_set")], top_abs_rank = direction_rank)[, c("target_name", "display_label", "branch_id", "model", "gene", "gene_order", "module", "predictor_class", "coefficient", "scaled_GBI", "contribution", "abs_contribution", "top_abs_rank", "top_set")]
)
top_out <- unique(top_out)
write_csv(top_out, file.path(out_dir, "SourceData_FigS3_ancestor_top_contributors.csv"))

module_rows <- do.call(rbind, lapply(split(fp, paste(fp$target_name, fp$model, fp$module, sep = "||")), function(part) {
  first <- part[1, , drop = FALSE]
  data.frame(
    target_name = first$target_name,
    display_label = first$display_label,
    branch_id = first$branch_id,
    branch_group = first$branch_group,
    model = first$model,
    module = first$module,
    net_contribution = sum(part$contribution),
    positive_contribution_sum = sum(part$contribution[part$contribution > 0]),
    negative_contribution_sum = sum(part$contribution[part$contribution < 0]),
    absolute_contribution_sum = sum(abs(part$contribution)),
    n_contributing_genes = nrow(part),
    n_positive_contribution_genes = sum(part$contribution > 0),
    n_negative_contribution_genes = sum(part$contribution < 0),
    top_genes_by_abs_contribution = paste(head(part$gene[order(part$abs_contribution, decreasing = TRUE)], 5), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))
write_csv(module_rows, file.path(out_dir, "SourceData_FigS3_ancestor_module_contributions.csv"))

message("[", Sys.time(), "] Building figures")
INK <- "#1F2430"
MUTED <- "#6F768A"
GRID <- "#E6E8F0"
MARINE_BLUE <- "#2C7FB8"
AQUATIC_ORANGE <- "#D95F02"
RED <- "#8B0000"
BLUE_NEG <- "#1F4E79"

score_plot_df <- profile_wide
score_plot_df$display_label <- factor(score_plot_df$display_label, levels = rev(profile_wide$display_label))
score_plot_df$label_x <- score_plot_df$marine_profile_score + 0.015
score_plot_df$label_y <- score_plot_df$aquatic_profile_score
manual_label_positions <- data.frame(
  display_label = c(
    "Crown Pinnipedia", "Otarioidea", "Otariidae", "Sirenia",
    "Dugong + Hydrodamalis", "Brown bear + polar bear",
    "Lutrinae / Enhydra-related"
  ),
  label_x = c(0.045, 0.070, 0.310, 0.050, 0.052, 0.052, 0.052),
  label_y = c(0.123, 0.147, 0.095, 0.205, 0.034, 0.002, 0.060),
  stringsAsFactors = FALSE
)
for (i in seq_len(nrow(manual_label_positions))) {
  idx <- as.character(score_plot_df$display_label) == manual_label_positions$display_label[i]
  score_plot_df$label_x[idx] <- manual_label_positions$label_x[i]
  score_plot_df$label_y[idx] <- manual_label_positions$label_y[i]
}
arrow_pairs <- data.frame(
  from = c(
    "Cetacea + hippo ancestor", "Crown Cetacea",
    "Crown Pinnipedia", "Crown Pinnipedia", "Otarioidea",
    "Sirenia"
  ),
  to = c(
    "Crown Cetacea", "Odontoceti",
    "Phocidae", "Otarioidea", "Otariidae",
    "Dugong + Hydrodamalis"
  ),
  stringsAsFactors = FALSE
)
arrow_df <- do.call(rbind, lapply(seq_len(nrow(arrow_pairs)), function(i) {
  from_row <- score_plot_df[as.character(score_plot_df$display_label) == arrow_pairs$from[i], , drop = FALSE]
  to_row <- score_plot_df[as.character(score_plot_df$display_label) == arrow_pairs$to[i], , drop = FALSE]
  if (nrow(from_row) == 0 || nrow(to_row) == 0) return(NULL)
  data.frame(
    x = from_row$marine_profile_score,
    y = from_row$aquatic_profile_score,
    xend = to_row$marine_profile_score,
    yend = to_row$aquatic_profile_score,
    stringsAsFactors = FALSE
  )
}))

plot_a <- ggplot(score_plot_df, aes(x = marine_profile_score, y = aquatic_profile_score)) +
  geom_segment(
    data = arrow_df,
    aes(x = x, y = y, xend = xend, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = unit(0.11, "inches")), color = "#B8BDC8", linewidth = 0.45
  ) +
  geom_point(aes(fill = branch_group), shape = 21, color = INK, size = 3.2, stroke = 0.35) +
  geom_text(aes(x = label_x, y = label_y, label = display_label), hjust = 0, vjust = 0.5, size = 2.35, color = INK) +
  scale_fill_manual(values = c(
    "Cetacean axis" = "#5DA5DA",
    "Pinniped axis" = "#9C755F",
    "Sirenian / edge context" = "#F1A340"
  )) +
  coord_cartesian(xlim = c(0, 1.12), ylim = c(0, 1.03), clip = "off") +
  labs(
    title = "A. Projected internal-branch profiles",
    subtitle = "Coordinates are descriptive profile scores under the Figure 6B/C full-data architectures",
    x = "Marine-axis profile score",
    y = "Aquatic-dependence profile score",
    fill = "Branch group"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    plot.title = element_text(face = "bold", color = INK, size = 11),
    plot.subtitle = element_text(color = MUTED, size = 7.5),
    axis.title = element_text(color = INK),
    axis.text = element_text(color = INK, size = 7.5),
    panel.grid = element_line(color = GRID, linewidth = 0.35),
    legend.position = "bottom",
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    plot.margin = margin(5.5, 30, 5.5, 5.5)
  )

label_top_genes <- function(d, n = 6) {
  do.call(rbind, lapply(split(d, paste(d$target_name, d$model, sep = "||")), function(part) {
    part <- part[order(part$abs_contribution, decreasing = TRUE), , drop = FALSE]
    part$label_rank <- seq_len(nrow(part))
    part$show_label <- part$label_rank <= n
    part
  }))
}
fp_labeled <- label_top_genes(fp, 6)

plot_fingerprints <- function(model_name, panel_title, n_expected) {
  d <- fp_labeled[fp_labeled$model == model_name, , drop = FALSE]
  d$display_label <- factor(d$display_label, levels = rev(targets$display_label))
  d$bar_color <- ifelse(d$contribution >= 0, "positive", "negative")
  ylim <- max(abs(d$contribution), na.rm = TRUE) * 1.24
  ggplot(d, aes(x = gene_order, y = contribution, fill = bar_color)) +
    geom_col(width = 0.82, color = NA, alpha = 0.94) +
    geom_hline(yintercept = 0, color = "#BFC5D0", linewidth = 0.35) +
    geom_text(
      data = subset(d, show_label),
      aes(label = gene, y = contribution + ifelse(contribution >= 0, ylim * 0.035, -ylim * 0.035)),
      angle = 90,
      size = 1.72,
      fontface = "bold",
      vjust = ifelse(subset(d, show_label)$contribution >= 0, 0, 1),
      color = ifelse(subset(d, show_label)$contribution >= 0, RED, BLUE_NEG)
    ) +
    facet_grid(display_label ~ ., switch = "y") +
    scale_fill_manual(values = c("positive" = RED, "negative" = BLUE_NEG), guide = "none") +
    coord_cartesian(ylim = c(-ylim, ylim), clip = "off") +
    labs(
      title = panel_title,
      subtitle = paste0("All ", n_expected, " nonzero predictors in fixed Figure 6 coefficient order; labels mark top absolute contributors"),
      x = "Genes in fixed full-data coefficient order",
      y = "Contribution to fitted logit"
    ) +
    theme_minimal(base_size = 8.5) +
    theme(
      plot.title = element_text(face = "bold", color = INK, size = 11),
      plot.subtitle = element_text(color = MUTED, size = 7.2),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(color = INK, size = 6.4),
      axis.title = element_text(color = INK, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "#EDF0F5", linewidth = 0.25),
      strip.text.y.left = element_text(angle = 0, hjust = 1, color = INK, face = "bold", size = 6.2),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.15, "lines"),
      plot.margin = margin(5.5, 9, 5.5, 14)
    )
}

plot_b <- plot_fingerprints("marine", "B. Marine-axis ancestral fingerprints", 71)
plot_c <- plot_fingerprints("binary_aquatic_dependence", "C. Aquatic-dependence ancestral fingerprints", 101)

module_heat <- module_rows
module_heat <- module_heat[module_heat$module != "Table S5-only / unassigned", , drop = FALSE]
top_modules <- names(sort(tapply(module_heat$absolute_contribution_sum, module_heat$module, sum), decreasing = TRUE))[1:min(12, length(unique(module_heat$module)))]
module_heat <- module_heat[module_heat$module %in% top_modules, , drop = FALSE]
module_heat$display_label <- factor(module_heat$display_label, levels = rev(targets$display_label))
module_heat$module <- factor(module_heat$module, levels = rev(top_modules))
module_heat$model_label <- ifelse(module_heat$model == "marine", "Marine axis", "Aquatic-dependence axis")
heat_lim <- max(abs(module_heat$net_contribution), na.rm = TRUE)
plot_d <- ggplot(module_heat, aes(x = module, y = display_label, fill = net_contribution)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_wrap(~ model_label, ncol = 1) +
  scale_fill_gradient2(low = BLUE_NEG, mid = "white", high = RED, midpoint = 0, limits = c(-heat_lim, heat_lim), name = "Net\ncontribution") +
  labs(
    title = "D. Functional-system contribution summaries",
    subtitle = "Descriptive module-annotated contribution sums; not module recurrence or enrichment",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 8.5) +
  theme(
    plot.title = element_text(face = "bold", color = INK, size = 11),
    plot.subtitle = element_text(color = MUTED, size = 7.2),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 5.8, color = INK),
    axis.text.y = element_text(size = 6.2, color = INK),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold", color = INK),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6)
  )

save_plot <- function(plot_obj, filename_base, width, height) {
  ggsave(file.path(out_dir, paste0(filename_base, ".pdf")), plot_obj, width = width, height = height, units = "in", device = grDevices::pdf)
  ggsave(file.path(out_dir, paste0(filename_base, ".png")), plot_obj, width = width, height = height, units = "in", dpi = 300)
  ggsave(file.path(out_dir, paste0(filename_base, ".svg")), plot_obj, width = width, height = height, units = "in", device = svglite::svglite)
}

pdf_file <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints.pdf")
grDevices::pdf(pdf_file, width = 15.6, height = 18.4, onefile = TRUE)
grid.newpage()
grid.draw(arrangeGrob(plot_a, plot_b, plot_c, plot_d, ncol = 2, heights = c(0.70, 1.30), top = textGrob("Supplementary Fig. S3 | Internal-branch profile fingerprints", gp = gpar(fontsize = 15, fontface = "bold", col = INK))))
dev.off()

png_file <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints.png")
png(png_file, width = 15.6, height = 18.4, units = "in", res = 300)
grid.newpage()
grid.draw(arrangeGrob(plot_a, plot_b, plot_c, plot_d, ncol = 2, heights = c(0.70, 1.30), top = textGrob("Supplementary Fig. S3 | Internal-branch profile fingerprints", gp = gpar(fontsize = 15, fontface = "bold", col = INK))))
dev.off()

svg_file <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints.svg")
svglite::svglite(svg_file, width = 15.6, height = 18.4)
grid.newpage()
grid.draw(arrangeGrob(plot_a, plot_b, plot_c, plot_d, ncol = 2, heights = c(0.70, 1.30), top = textGrob("Supplementary Fig. S3 | Internal-branch profile fingerprints", gp = gpar(fontsize = 15, fontface = "bold", col = INK))))
dev.off()

multi_pdf <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_multipage_all_nodes.pdf")
grDevices::pdf(multi_pdf, width = 11.8, height = 14.0, onefile = TRUE)
print(plot_a)
print(plot_b)
print(plot_c)
print(plot_d)
dev.off()

message("[", Sys.time(), "] Writing QC and narrative notes")
qc_lines <- c(
  "# Ancestor Fingerprint QC Report",
  "",
  "## Evidence Layer",
  "This supplement uses the Figure 6B/C corrected full-data architecture: 71 marine predictors and 101 binary aquatic-dependence predictors. Contributions equal scaled GBI multiplied by the corrected full-data coefficient; intercepts are omitted from barplots but included in profile-score reconstruction.",
  "",
  "Internal branches are descriptive genomic-profile projection targets, not direct ancestral habitat assignments and not out-of-fold validation.",
  "",
  "## Full-Data Model Rebuild",
  paste(capture.output(print(coef_repro, row.names = FALSE)), collapse = "\n"),
  "",
  "## Branch Target Check",
  paste(capture.output(print(targets[, c("target_name", "display_label", "branch_id", "target_type", "inclusion_set", "branch_id_found_in_gbi")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Reconstructed Profile Scores",
  paste(capture.output(print(profile_wide[, c("display_label", "branch_id", "marine_profile_score", "aquatic_profile_score", "marine_n_predictors", "aquatic_n_predictors", "figure6A_source_comparison_status")], row.names = FALSE)), collapse = "\n"),
  "",
  "## Figure 6A Source Comparison",
  "Current Figure 6A internal rows use terminal-only internal branch projections, whereas this supplement uses the Figure 6B/C corrected full-data 71/101 architecture. Differences from the current Figure 6A source are therefore recorded as evidence-layer differences rather than treated as numerical errors.",
  "",
  paste0("- Maximum absolute logit reconstruction error across model/branch rows: ", format(max(profile_rows$abs_logit_reconstruction_error), scientific = TRUE)),
  paste0("- Marine predictor count used in fingerprints: ", length(models$marine$selected)),
  paste0("- Aquatic predictor count used in fingerprints: ", length(models$binary_aquatic_dependence$selected))
)
writeLines(qc_lines, file.path(out_dir, "ancestor_fingerprint_QC_report.md"))

caption <- c(
  "**Supplementary Fig. S3 | Gene-level fingerprints underlying projected ancestral-branch profiles.**",
  "(A) Selected internal branches projected onto the marine and binary aquatic-dependence full-data profile axes. Coordinates are descriptive genomic-profile scores under the same corrected full-data architectures used for Figure 6B/C and should not be interpreted as direct ancestral habitat assignments.",
  "(B,C) Marine-axis (B) and binary aquatic-dependence (C) fingerprints for selected internal branches. Bars show all nonzero predictors in the fixed Figure 6 coefficient order, with contributions calculated as scaled GBI x corrected full-data LASSO coefficient; model intercepts are omitted from bars but included in profile-score reconstruction. Gene labels mark the largest absolute contributors in each branch.",
  "(D) Descriptive module-annotated contribution summaries aggregate predictor contributions by Figure 4C/Supplementary Table S5 display module. These summaries are not pathway enrichment tests and do not establish module-level recurrence or convergence."
)
writeLines(caption, file.path(out_dir, "ancestor_fingerprint_caption_draft.md"))

top_summary <- function(target, model_name, n = 6) {
  part <- fp[fp$target_name == target & fp$model == model_name, , drop = FALSE]
  part <- part[order(part$abs_contribution, decreasing = TRUE), , drop = FALSE]
  paste(head(part$gene, n), collapse = ", ")
}
module_summary <- function(target, model_name, n = 3) {
  part <- module_rows[module_rows$target_name == target & module_rows$model == model_name, , drop = FALSE]
  part <- part[order(part$absolute_contribution_sum, decreasing = TRUE), , drop = FALSE]
  paste(head(part$module, n), collapse = "; ")
}
cet_delta <- profile_wide$marine_profile_score[profile_wide$target_name == "crown_Cetacea"] - profile_wide$marine_profile_score[profile_wide$target_name == "common_ancestor_Cetacea_plus_Hippopotamus"]
odont_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Odontoceti"]
odont_a <- profile_wide$aquatic_profile_score[profile_wide$target_name == "Odontoceti"]
pinn_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Pinnipedia"]
phoc_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Phocidae"]
otar_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Otarioidea"]
siren_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Sirenia"]
siren_a <- profile_wide$aquatic_profile_score[profile_wide$target_name == "Sirenia"]
dug_m <- profile_wide$marine_profile_score[profile_wide$target_name == "Dugong_plus_Hydrodamalis"]
dug_a <- profile_wide$aquatic_profile_score[profile_wide$target_name == "Dugong_plus_Hydrodamalis"]

story <- c(
  "# Section 2.6 Ancestor Gene Story Suggestions",
  "",
  "## Computed from Figure 6B/C 71/101 corrected full-data architecture",
  paste0("- Crown Cetacea marine score minus Cetacea + hippo ancestor marine score: ", fmt3(cet_delta), "."),
  paste0("- Crown Cetacea top marine-axis genes: ", top_summary("crown_Cetacea", "marine"), "."),
  paste0("- Crown Cetacea top marine-axis modules: ", module_summary("crown_Cetacea", "marine"), "."),
  paste0("- Odontoceti profile scores: marine ", fmt3(odont_m), ", aquatic ", fmt3(odont_a), "."),
  paste0("- Odontoceti top marine-axis genes: ", top_summary("Odontoceti", "marine"), "."),
  paste0("- Odontoceti top aquatic-axis genes: ", top_summary("Odontoceti", "binary_aquatic_dependence"), "."),
  paste0("- Crown Pinnipedia marine score: ", fmt3(pinn_m), "; Phocidae: ", fmt3(phoc_m), "; Otarioidea: ", fmt3(otar_m), "."),
  paste0("- Phocidae top marine-axis genes: ", top_summary("Phocidae", "marine"), "."),
  paste0("- Otarioidea top marine-axis genes: ", top_summary("Otarioidea", "marine"), "."),
  paste0("- Sirenia profile scores: marine ", fmt3(siren_m), ", aquatic ", fmt3(siren_a), "; Dugong + Hydrodamalis: marine ", fmt3(dug_m), ", aquatic ", fmt3(dug_a), "."),
  paste0("- Sirenia top aquatic-axis genes: ", top_summary("Sirenia", "binary_aquatic_dependence"), "."),
  "",
  "## Insertable Section 2.6 sentences",
  paste0("Projected internal branches showed that the shift from the Cetacea + hippo ancestor to crown Cetacea involved a broader change in fitted genomic-profile coordinates, with the crown Cetacea fingerprint driven by contributors such as ", top_summary("crown_Cetacea", "marine", 4), " on the marine axis."),
  paste0("Odontoceti retained a distinct projected fingerprint on both axes (marine score ", fmt3(odont_m), "; aquatic-dependence score ", fmt3(odont_a), "), with top contributions distributed across gene-level components rather than a single predictor."),
  paste0("In contrast, crown Pinnipedia had a weak fitted projection in this internal-branch analysis (marine score ", fmt3(pinn_m), "), whereas Phocidae showed a stronger marine-axis fingerprint than Otarioidea, consistent with lineage-weighted routes within pinnipeds."),
  paste0("Sirenia and the Dugong + Hydrodamalis branch showed partial or distinct profiles (Sirenia marine/aquatic scores ", fmt3(siren_m), "/", fmt3(siren_a), "; Dugong + Hydrodamalis ", fmt3(dug_m), "/", fmt3(dug_a), "), supporting a decoupling interpretation rather than a uniform marine-like route."),
  "These branch-level profiles are descriptive genomic projections built from scaled GBI x coefficient contributions; they are not ancestral habitat reconstructions or gene-level functional validation."
)
writeLines(story, file.path(out_dir, "Section2_6_ancestor_gene_story_suggestions.md"))

readme <- c(
  "# Ancestor Fingerprints Supplement",
  "",
  "This directory was generated by `code/12_NC_species_ancestral_fingerprints/build_FigureS3_ancestor_fingerprints.R`.",
  "",
  "No GBI, ASR, t-test/FDR discovery, enrichment, or LASSO feature-selection analyses were rerun. The script rebuilds the two corrected full-data baseline LASSO fits only to recover intercepts and full-data imputation/scaling parameters, then verifies exact coefficient reproduction against Phase 12A archived coefficients before projecting selected internal branches.",
  "",
  "Main outputs:",
  "- `SupplementaryFig_ancestor_fingerprints.pdf/png/svg`",
  "- `SupplementaryFig_ancestor_fingerprints_multipage_all_nodes.pdf`",
  "- `SourceData_FigS3_ancestor_profile_scores.csv`",
  "- `SourceData_FigS3_ancestor_fingerprints_long.csv`",
  "- `SourceData_FigS3_ancestor_top_contributors.csv`",
  "- `SourceData_FigS3_ancestor_module_contributions.csv`",
  "- `ancestor_fingerprint_QC_report.md`",
  "- `ancestor_fingerprint_caption_draft.md`",
  "- `Section2_6_ancestor_gene_story_suggestions.md`",
  "",
  "Evidence boundary: internal branches are descriptive projected genomic profiles, not direct ancestral habitat assignments."
)
writeLines(readme, file.path(out_dir, "README_folded_ancestor_fingerprints_supplement.md"))

manifest_files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
manifest_files <- manifest_files[file.info(manifest_files)$isdir == FALSE]
sha <- vapply(manifest_files, function(f) system2("shasum", c("-a", "256", f), stdout = TRUE)[1], character(1))
manifest <- data.frame(
  file = sub(paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", out_dir), "/?"), "", manifest_files),
  size_bytes = file.info(manifest_files)$size,
  sha256 = sub(" .*", "", sha),
  stringsAsFactors = FALSE
)
write_tsv(manifest[order(manifest$file), ], file.path(out_dir, "file_manifest_with_sha256.tsv"))

message("[", Sys.time(), "] Done: ", out_dir)
