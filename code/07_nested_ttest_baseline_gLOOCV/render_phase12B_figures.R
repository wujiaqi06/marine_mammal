#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(svglite)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "10_reviewer_risk_controls", "12B_nested_ttest_baseline_robustness")
out_dir <- file.path(work_dir, "nested_ttest_baseline")
fig_dir <- file.path(work_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(file, ...) {
  utils::read.delim(file, stringsAsFactors = FALSE, check.names = FALSE, ...)
}

save_plot_all <- function(plot, stem, width, height) {
  grDevices::pdf(file.path(fig_dir, paste0(stem, ".pdf")), width = width, height = height, useDingbats = FALSE)
  print(plot)
  grDevices::dev.off()
  grDevices::png(file.path(fig_dir, paste0(stem, ".png")), width = width, height = height, units = "in", res = 300)
  print(plot)
  grDevices::dev.off()
  svglite::svglite(file.path(fig_dir, paste0(stem, ".svg")), width = width, height = height)
  print(plot)
  grDevices::dev.off()
}

roc_all <- read_tsv(file.path(out_dir, "nested_ttest_ROC_points.tsv"))
auc_summary <- read_tsv(file.path(out_dir, "nested_ttest_AUC_summary.tsv"))
oof <- read_tsv(file.path(out_dir, "nested_ttest_OOF_predictions.tsv"))

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
  theme(
    legend.position = c(0.62, 0.18),
    legend.background = element_rect(fill = "white", color = "grey85"),
    plot.title = element_text(face = "bold", size = 10)
  )
save_plot_all(p_roc, "Figure2B_nested_ttest_ROC_candidate", 4.2, 3.4)

dist_df <- oof[oof$nested_OOF_available == "TRUE" | oof$nested_OOF_available == TRUE, , drop = FALSE]
dist_df$trait <- ifelse(grepl("marine_binary", dist_df$model), "marine_binary", "binary_aquatic_dependence")
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
      x = 2.0,
      y = 0.085,
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
svglite::svglite(file.path(fig_dir, "Figure2_BC_nested_ttest_candidate_combined.svg"), width = 10.5, height = 3.7)
grid::grid.draw(p_combined)
grDevices::dev.off()

message("Phase 12B figures rendered with svglite SVG output.")
