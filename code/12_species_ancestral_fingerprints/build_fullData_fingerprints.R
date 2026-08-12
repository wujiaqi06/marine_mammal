#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glmnet)
})

root_value <- Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "")
if (!nzchar(root_value)) {
  stop("Set MARINE_MAMMAL_ENDPOINTFIX_ROOT to the unpacked endpoint-fix analysis root.")
}
root_dir <- normalizePath(root_value, mustWork = TRUE)
output_root <- Sys.getenv(
  "MARINE_MAMMAL_OUTPUT_ROOT",
  unset = file.path(root_dir, "reproduction_outputs")
)
archive_root <- Sys.getenv("MARINE_MAMMAL_ARCHIVE_ROOT", unset = output_root)
out_dir <- file.path(output_root, "figure6_fullData_fingerprints")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(file, ...) utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
write_csv <- function(x, file) utils::write.csv(x, file = file, quote = TRUE, row.names = FALSE, na = "")
write_tsv <- function(x, file) utils::write.table(x, file = file, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
inv_logit <- function(x) 1 / (1 + exp(-x))

paths <- list(
  gbi = file.path(root_dir, "04_GBI_matrix", "branch_label_crosswalk", "outputs", "endpointfix_no_fuse.fix.GBI_matrix.oldlabels.tsv"),
  branch = file.path(root_dir, "07_ttest_screening", "inputs", "branch_files", "mammal.branch.txt"),
  trait = file.path(root_dir, "05_trait_tables", "derived", "trait_table.mammal302.active_TY_NK_final_18pt.pipeline_alias.tsv"),
  marine_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_marine_binary", "marine_binary.mammal.FDR0.01.n1559.t_test.txt"),
  aquatic_features = file.path(root_dir, "07_ttest_screening", "stage1_deterministic_asr", "outputs", "fix_aquatic_v2", "aquatic_v2.mammal.FDR0.01.n1227.t_test.txt"),
  archived_coef = file.path(root_dir, "10_reviewer_risk_controls", "12A_corrected_preprocessing_LASSO_architecture_sensitivity", "corrected_coefficients", "corrected_preprocessing_full_data_coefficients_all_runs.tsv"),
  fig4_architecture = file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment", "Figure4_corrected_full_data", "Figure4_corrected_coefficient_architecture.tsv"),
  table_s5 = file.path(root_dir, "supplementary_tables_endpointfix", "TableS5", "Supplementary_Table_S5_Figure4C_predictor_annotation.tsv"),
  terminal_oof = file.path(archive_root, "terminal_oof_predictions.csv")
)
missing <- names(paths)[!file.exists(unlist(paths))]
if (length(missing) > 0) stop("Missing required inputs: ", paste(missing, collapse = ", "))

message("[", Sys.time(), "] Reading frozen inputs")
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
fig4_arch <- read_tsv(paths$fig4_architecture)
table_s5 <- read_tsv(paths$table_s5)
terminal_oof <- utils::read.csv(paths$terminal_oof, stringsAsFactors = FALSE, check.names = FALSE)

display_label <- function(species) {
  parts <- strsplit(species, "_", fixed = TRUE)[[1]]
  if (length(parts) >= 2) paste(parts[1], parts[2]) else species
}

clade_map <- c(
  Orcinus_orca = "cetacean",
  Zalophus_californianus = "pinniped",
  Leptonychotes_weddellii = "pinniped",
  Odobenus_rosmarus_divergens = "pinniped",
  Ursus_maritimus = "marine edge",
  Enhydra_lutris_kenyoni = "marine edge",
  Dugong_dugon = "sirenian",
  Trichechus_manatus_latirostris = "sirenian",
  Hydrodamalis_gigas = "sirenian",
  Hippopotamus_amphibius = "non-marine aquatic control",
  Lutra_lutra = "non-marine aquatic control",
  Lontra_canadensis = "non-marine aquatic control",
  Pteronura_brasiliensis = "non-marine aquatic control",
  Aonyx_cinereus = "non-marine aquatic control",
  Platanista_minor = "freshwater cetacean",
  Inia_geoffrensis = "freshwater cetacean",
  Lipotes_vexillifer = "freshwater cetacean"
)
focal_species <- names(clade_map)
missing_focal <- setdiff(focal_species, terminal$Species)
if (length(missing_focal) > 0) warning("Missing focal species in terminal table: ", paste(missing_focal, collapse = ", "))
focal_species <- intersect(focal_species, terminal$Species)

make_spec <- function(model) {
  if (model == "marine") {
    genes <- read_tsv(paths$marine_features)$gene
    list(
      model = "marine",
      run_id = "fix_marine_binary",
      genes = genes,
      y = as.numeric(terminal$marine_binary),
      eligible = is.finite(terminal$marine_binary),
      n_expected_selected = 71L
    )
  } else {
    genes <- read_tsv(paths$aquatic_features)$gene
    y <- as.numeric(terminal$aquatic_v2)
    list(
      model = "binary_aquatic_dependence",
      run_id = "fix_aquatic_v2",
      genes = genes,
      y = y,
      eligible = is.finite(y) & y != 0.5,
      n_expected_selected = 101L
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
  if (ncol(x_raw) == 0) stop("No usable features for ", spec$model)
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
  intercept <- unname(coef_all["(Intercept)"])
  selected <- names(beta)[beta != 0]
  list(
    spec = spec,
    eligible_species = eligible_species,
    x_scaled = x_scaled,
    y = y,
    impute_means = impute_means,
    scale_means = scale_means,
    scale_sds = scale_sds,
    lambda_min = cv$lambda.min,
    intercept = intercept,
    beta = beta,
    selected = selected,
    dropped_all_missing = sum(!keep_not_all_missing),
    dropped_zero_variance = sum(!keep_nonzero)
  )
}

project_species <- function(model_fit, species) {
  spec <- model_fit$spec
  row <- terminal[match(species, terminal$Species), , drop = FALSE]
  genes <- model_fit$selected
  raw <- as.numeric(gbi_mat[genes, row$Branch])
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
    species = species,
    display_label = display_label(species),
    clade_group = unname(clade_map[species]),
    trait_category = row$trait_category,
    aquaticity_score = row$aquaticity_score,
    model = spec$model,
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
    stringsAsFactors = FALSE
  )
}

add_annotation <- function(df) {
  module_col <- if ("display_module" %in% names(table_s5)) "display_module" else if ("candidate_module" %in% names(table_s5)) "candidate_module" else NA_character_
  predictor_cols <- c("gene", "lasso_group", "display_module", "pro_recommended_module", "recommended_submodule_or_function", "annotation_confidence")
  predictor_cols <- intersect(predictor_cols, names(table_s5))
  ann <- table_s5[, predictor_cols, drop = FALSE]
  names(ann)[names(ann) == "lasso_group"] <- "predictor_class"
  if (!"display_module" %in% names(ann)) ann$display_module <- NA_character_
  if (!"predictor_class" %in% names(ann)) ann$predictor_class <- NA_character_
  merge(df, ann, by = "gene", all.x = TRUE, sort = FALSE)
}

message("[", Sys.time(), "] Rebuilding corrected full-data baseline LASSO fits")
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
write_tsv(coef_repro, file.path(out_dir, "fullData_model_rebuild_coefficient_reproducibility.tsv"))
if (any(coef_repro$status != "PASS")) {
  stop("Full-data model rebuild did not reproduce archived selected coefficients; see fullData_model_rebuild_coefficient_reproducibility.tsv")
}

message("[", Sys.time(), "] Projecting focal species into corrected full-data coordinate systems")
fp <- do.call(rbind, c(
  lapply(focal_species, function(s) project_species(models$marine, s)),
  lapply(focal_species, function(s) project_species(models$binary_aquatic_dependence, s))
))
fp <- add_annotation(fp)
fp$predictor_class <- ifelse(is.na(fp$predictor_class) | fp$predictor_class == "", "not_classified", fp$predictor_class)
fp$module <- ifelse(!is.na(fp$display_module) & fp$display_module != "", fp$display_module, "Table S5-only / unassigned")

gene_order_list <- lapply(models, function(m) {
  genes <- m$selected
  coefs <- m$beta[genes]
  ord <- order(coefs, decreasing = FALSE)
  genes <- genes[ord]
  tmp <- data.frame(
    model = m$spec$model,
    gene = genes,
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
fp <- fp[order(match(fp$model, c("marine", "binary_aquatic_dependence")), match(fp$species, focal_species), fp$gene_order), , drop = FALSE]
fingerprint_out <- fp[, c(
  "species", "display_label", "clade_group", "trait_category", "aquaticity_score",
  "model", "gene", "gene_order", "predictor_class", "module", "coefficient",
  "scaled_GBI", "contribution", "abs_contribution", "contribution_sign",
  "GBI_raw_if_available", "imputed_flag"
), drop = FALSE]
write_csv(fingerprint_out, file.path(out_dir, "SourceData_Fig6BC_species_fingerprints_long.csv"))
write_csv(gene_order_out, file.path(out_dir, "Fig6_fullData_gene_order_map.csv"))

projection_rows <- do.call(rbind, lapply(split(fp, paste(fp$species, fp$model, sep = "||")), function(part) {
  first <- part[1, , drop = FALSE]
  oof_prob <- NA_real_
  oof_logit <- NA_real_
  if (first$model == "marine" && first$species %in% terminal_oof$species) {
    oof_prob <- terminal_oof$marine_OOF_probability[match(first$species, terminal_oof$species)]
    oof_logit <- terminal_oof$marine_OOF_logit[match(first$species, terminal_oof$species)]
  }
  if (first$model == "binary_aquatic_dependence" && first$species %in% terminal_oof$species) {
    oof_prob <- terminal_oof$binary_aquatic_dependence_OOF_probability[match(first$species, terminal_oof$species)]
    oof_logit <- terminal_oof$binary_aquatic_dependence_OOF_logit[match(first$species, terminal_oof$species)]
  }
  data.frame(
    species = first$species,
    display_label = first$display_label,
    clade_group = first$clade_group,
    trait_category = first$trait_category,
    aquaticity_score = first$aquaticity_score,
    model = first$model,
    n_predictors = nrow(part),
    n_imputed_predictors = sum(part$imputed_flag),
    intercept = models[[first$model]]$intercept,
    contribution_sum = sum(part$contribution),
    fitted_logit = first$fitted_logit,
    fitted_probability = first$fitted_probability,
    OOF_probability_if_available = oof_prob,
    OOF_logit_if_available = oof_logit,
    fitted_minus_OOF = first$fitted_probability - oof_prob,
    evidence_label = "descriptive full-data fitted projection after corrected preprocessing; not held-out validation",
    stringsAsFactors = FALSE
  )
}))
projection_rows <- projection_rows[order(match(projection_rows$species, focal_species), projection_rows$model), , drop = FALSE]
write_csv(projection_rows, file.path(out_dir, "SourceData_Fig6_fullData_species_projections_all_focal.csv"))

logit_check <- do.call(rbind, lapply(split(fp, paste(fp$species, fp$model, sep = "||")), function(part) {
  first <- part[1, , drop = FALSE]
  model <- models[[first$model]]
  reconstructed <- model$intercept + sum(part$contribution)
  data.frame(
    species = first$species,
    model = first$model,
    n_predictors = nrow(part),
    intercept = model$intercept,
    sum_contribution = sum(part$contribution),
    stored_fitted_logit = first$fitted_logit,
    reconstructed_logit = reconstructed,
    abs_logit_difference = abs(first$fitted_logit - reconstructed),
    status = ifelse(abs(first$fitted_logit - reconstructed) < 1e-12, "PASS", "FAIL"),
    stringsAsFactors = FALSE
  )
}))
write_tsv(logit_check, file.path(out_dir, "fullData_logit_reconstruction_check.tsv"))

impute_summary <- aggregate(
  imputed_flag ~ species + display_label + clade_group + model,
  fingerprint_out,
  function(x) sum(x)
)
names(impute_summary)[names(impute_summary) == "imputed_flag"] <- "n_imputed_predictors"
impute_summary$n_predictors <- ifelse(impute_summary$model == "marine", length(models$marine$selected), length(models$binary_aquatic_dependence$selected))
write_tsv(impute_summary, file.path(out_dir, "fullData_imputation_count_summary.tsv"))

writeLines(c(
  "# Figure 6 Full-Data QC Report",
  "",
  "## Located Input Files",
  paste0("- GBI matrix: `", paths$gbi, "`"),
  paste0("- Branch table: `", paths$branch, "`"),
  paste0("- Trait table: `", paths$trait, "`"),
  paste0("- Marine candidate genes: `", paths$marine_features, "`"),
  paste0("- Binary aquatic candidate genes: `", paths$aquatic_features, "`"),
  paste0("- Archived Phase 12A coefficients: `", paths$archived_coef, "`"),
  paste0("- Supplementary Table S5: `", paths$table_s5, "`"),
  paste0("- OOF comparison table: `", paths$terminal_oof, "`"),
  "",
  "## Model Rebuild",
  paste(capture.output(print(coef_repro, row.names = FALSE)), collapse = "\n"),
  "",
  "The full-data models were rebuilt only to recover intercepts and full-data imputation/scaling parameters that were not archived in the coefficient table. The rebuilt selected gene sets and coefficients exactly match the archived Phase 12A corrected full-data baseline coefficients within numerical tolerance.",
  "",
  "## Predictor Counts",
  paste0("- Marine full-data predictors: ", length(models$marine$selected)),
  paste0("- Binary aquatic-dependence full-data predictors: ", length(models$binary_aquatic_dependence$selected)),
  "",
  "## Projection Boundary",
  "Full-data fitted projections are descriptive architecture summaries after corrected preprocessing. They are not held-out validation estimates. Nested gLOOCV validation remains Fig. 2/Fig. 5A; fold-specific held-out decompositions are separate sensitivity material.",
  "",
  "## Logit Reconstruction",
  paste0("- Maximum absolute logit reconstruction error: ", format(max(logit_check$abs_logit_difference), scientific = TRUE)),
  "",
  "## Imputation Summary",
  paste(capture.output(print(impute_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "## Full-Data vs OOF Notes",
  "The species projection table reports OOF probabilities when available and fitted-minus-OOF differences for claim-boundary review. Large differences should be framed as fitted architecture/OOF discordance, not as stronger species-level validation."
), file.path(out_dir, "Fig6_fullData_QC_report.md"))

claim_lines <- c(
  "# Figure 6 Full-Data Claim Ceiling Report",
  "",
  "## Evidence Type",
  "Figure 6 full-data fingerprints use a shared corrected full-data LASSO coordinate system: 71 marine predictors and 101 binary aquatic-dependence predictors. These panels summarize fitted architecture after validation and should not be described as held-out validation.",
  "",
  "## Safe Main-Text Portraits",
  "- Orcinus orca: strong canonical cetacean portrait; use as a high marine-like and high aquatic-like fitted projection, cross-checked against OOF values.",
  "- Zalophus californianus: strong pinniped portrait suitable for main text.",
  "",
  "## Moderate / Heterogeneous Portraits",
  "- Leptonychotes weddellii and Odobenus rosmarus divergens should be framed as pinniped heterogeneity rather than uniform pinniped signal.",
  "",
  "## Edge / Decoupling Cases",
  "- Ursus maritimus, Enhydra lutris kenyoni and sirenians should be framed as ecological marine membership or aquatic dependence that can decouple from sparse marine-like fitted profiles.",
  "",
  "## River Dolphin Bridge",
  "- Platanista minor, Inia geoffrensis and Lipotes vexillifer are useful river-dolphin bridge cases for a fuller or extended version.",
  "",
  "## Avoid",
  "- Do not use posterior probability.",
  "- Do not call fitted full-data projections held-out validation.",
  "- Do not claim gene-level functional validation from contribution bars.",
  "- Do not describe edge taxa as confident marine predictions if fitted/OOF evidence is weak or discordant.",
  "",
  "## Recommended Caption Boundary",
  "Fitted projections are descriptive full-data architecture summaries after validation, not held-out performance estimates. Cross-validated performance is shown in Fig. 2 and Fig. 5A; fold-specific held-out decompositions are provided as Extended Data."
)
writeLines(claim_lines, file.path(out_dir, "Fig6_fullData_claim_ceiling_report.md"))

writeLines(c(
  "# Figure 6 Full-Data Fingerprints",
  "",
  "This directory contains Figure 6 full-data fitted projection drafts and source data.",
  "",
  "No GBI, ASR, t-test/FDR discovery, enrichment or nested validation analyses were rerun. The script rebuilds only the two corrected full-data baseline LASSO fits to recover intercepts and preprocessing parameters missing from archived coefficient outputs, and verifies exact coefficient reproduction against Phase 12A before exporting projections.",
  "",
  "Main evidence boundary: Fig. 4 full-data projections are descriptive fitted architecture summaries after validation. They are not held-out validation estimates."
), file.path(out_dir, "README_Fig6_fullData_fingerprints.md"))

message("[", Sys.time(), "] Wrote full-data fingerprint source data to ", out_dir)
