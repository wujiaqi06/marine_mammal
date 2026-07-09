#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(svglite)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
base_dir <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment")
fig_dir <- file.path(base_dir, "figures_draft")

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

save_svg <- function(plot, stem, width, height) {
  svglite::svglite(file.path(fig_dir, paste0(stem, ".svg")), width = width, height = height)
  print(plot)
  grDevices::dev.off()
}

theme_phase13 <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom"
    )
}

partition <- read_tsv(file.path(base_dir, "Figure4_corrected_full_data", "Figure4_corrected_predictor_partition.tsv"))
partition_plot <- subset(partition, partition %in% c("marine-specific predictors", "shared predictors", "aquatic-specific predictors"))
partition_plot$partition <- factor(
  partition_plot$partition,
  levels = c("marine-specific predictors", "shared predictors", "aquatic-specific predictors")
)
p4a <- ggplot(partition_plot, aes(partition, n_genes, fill = partition)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = n_genes), vjust = -0.35, size = 3.2) +
  scale_fill_manual(values = c("#2C7FB8", "#6A6A6A", "#F28E2B")) +
  labs(title = "Corrected full-data predictor partition", x = NULL, y = "Predictors") +
  theme_phase13() +
  theme(legend.position = "none")
save_svg(p4a, "Figure4A_corrected_full_data_partition_candidate", 4.8, 3.4)

coef_tab <- read_tsv(file.path(base_dir, "Figure4_corrected_full_data", "Figure4_corrected_coefficient_architecture.tsv"))
coef_top <- coef_tab[order(-coef_tab$max_abs_coef), ]
coef_top <- head(coef_top, 45)
coef_long <- rbind(
  data.frame(gene = coef_top$gene, model = "Marine", coefficient = coef_top$coefficient_marine, predictor_class = coef_top$display_category),
  data.frame(gene = coef_top$gene, model = "Aquatic", coefficient = coef_top$coefficient_aquatic, predictor_class = coef_top$display_category)
)
coef_long$gene <- factor(coef_long$gene, levels = rev(unique(coef_top$gene)))
p4b <- ggplot(coef_long, aes(gene, coefficient, fill = model)) +
  geom_col(position = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c(Marine = "#2C7FB8", Aquatic = "#F28E2B")) +
  labs(title = "Largest corrected full-data coefficients", x = NULL, y = "Scaled coefficient") +
  theme_phase13(7)
save_svg(p4b, "Figure4B_corrected_full_data_coefficients_candidate", 5.2, 3.6)

module_tab <- read_tsv(file.path(base_dir, "Figure4_corrected_full_data", "Figure4_corrected_module_partition.tsv"))
p4c <- ggplot(module_tab, aes(candidate_module, gene_count, fill = display_category)) +
  geom_col(position = "stack", color = "grey25", linewidth = 0.15) +
  coord_flip() +
  scale_fill_manual(values = c(
    "Marine-specific predictors" = "#2C7FB8",
    "Shared predictors" = "#6A6A6A",
    "Aquatic-specific predictors" = "#F28E2B"
  )) +
  labs(title = "Corrected predictor module partition", x = NULL, y = "Genes") +
  theme_phase13(7)
save_svg(p4c, "Figure4C_corrected_full_data_modules_candidate", 7.4, 4.2)

svglite::svglite(file.path(fig_dir, "Figure4_corrected_full_data_candidate.svg"), width = 8.6, height = 7.2)
gridExtra::grid.arrange(p4a, p4b, p4c, ncol = 1, heights = c(1.0, 1.25, 1.25))
grDevices::dev.off()

fig5a <- read_tsv(file.path(base_dir, "Figure5A_nested_sensitivity", "nested_ttest_Fig5A_sensitivity_AUC_summary.tsv"))
fig5a$run_id <- factor(fig5a$run_id, levels = rev(fig5a$run_id))
p5a <- ggplot(fig5a, aes(run_id, AUC_nested_ttest, fill = trait_axis)) +
  geom_col(width = 0.72, color = "grey20", linewidth = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey40") +
  geom_text(aes(label = sprintf("%.3f", AUC_nested_ttest)), hjust = -0.08, size = 2.7) +
  coord_flip(ylim = c(0, 1.05)) +
  scale_fill_manual(values = c("marine" = "#2C7FB8", "binary aquatic-dependence" = "#F28E2B")) +
  labs(title = "Nested t-test LASSO sensitivity", x = NULL, y = "AUC") +
  theme_phase13(8)
save_svg(p5a, "Figure5A_nested_ttest_candidate", 7.5, 4.8)

fig5b <- read_tsv(file.path(base_dir, "Figure5B_permutation_cleanup", "Figure5B_endpointfix_source_check.tsv"))
fig5b_plot <- data.frame(
  metric = c("Observed FDR genes", "Null median FDR genes"),
  value = c(fig5b$observed_n_sig_FDR_0_01[1], fig5b$null_median_sig[1])
)
p5b <- ggplot(fig5b_plot, aes(metric, value, fill = metric)) +
  geom_col(width = 0.62, color = "grey20", linewidth = 0.2) +
  geom_text(aes(label = value), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c("Observed FDR genes" = "#2C7FB8", "Null median FDR genes" = "#B9B9B9")) +
  labs(title = "Exact positive-count permutation control", x = NULL, y = "FDR-significant genes") +
  theme_phase13() +
  theme(legend.position = "none")
save_svg(p5b, "Figure5B_endpointfix_cleaned_candidate", 5.2, 3.6)

fig5c <- read_tsv(file.path(base_dir, "Figure5C_corrected_turnover_module_null", "Figure5C_module_null_summary.tsv"))
fig5c$comparison <- factor(fig5c$comparison, levels = unique(fig5c$comparison))
p5c <- ggplot(fig5c, aes(metric, observed_value, color = comparison)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.45)) +
  geom_point(aes(y = null_median), shape = 1, size = 2.2, position = position_dodge(width = 0.45)) +
  facet_wrap(~ comparison, ncol = 1) +
  labs(title = "Corrected turnover and candidate-union module null", x = NULL, y = "Metric value") +
  theme_phase13(7) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 25, hjust = 1))
save_svg(p5c, "Figure5C_corrected_turnover_with_module_null_candidate", 7.2, 4.2)

cat("SVG outputs rendered with svglite\n")
