#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(svglite)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
base_dir <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment")
fig_dir <- file.path(base_dir, "figures_draft")
fig4_dir <- file.path(base_dir, "Figure4_corrected_full_data")

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)

save_combined <- function(plot_list, stem, width = 7.6, height = 9.2) {
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  svg_path <- file.path(fig_dir, paste0(stem, ".svg"))

  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  grDevices::dev.off()

  grDevices::png(png_path, width = width, height = height, units = "in", res = 300)
  gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  grDevices::dev.off()

  svglite::svglite(svg_path, width = width, height = height)
  gridExtra::grid.arrange(grobs = plot_list, ncol = 1)
  grDevices::dev.off()
}

save_panel <- function(plot, stem, width, height) {
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

partition <- read_tsv(file.path(fig4_dir, "Figure4_corrected_predictor_partition.tsv"))
coef_arch <- read_tsv(file.path(fig4_dir, "Figure4_corrected_coefficient_architecture.tsv"))
module_partition <- read_tsv(file.path(fig4_dir, "Figure4_corrected_module_partition.tsv"))

p4a <- ggplot(partition[partition$partition != "union total", ], aes(x = partition, y = n_genes, fill = partition)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = n_genes), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("marine-specific predictors" = "#2F72C4", "shared predictors" = "#666666", "aquatic-specific predictors" = "#E28416")) +
  labs(x = NULL, y = "Predictors", title = "A. Corrected full-data predictor partition") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

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
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

p4c_df <- module_partition[module_partition$candidate_module != "REVIEW_UNANNOTATED", , drop = FALSE]
p4c <- ggplot(p4c_df, aes(x = display_category, y = gene_count, fill = candidate_module)) +
  geom_col(width = 0.72) +
  labs(x = NULL, y = "Module-annotated predictors", fill = "Module", title = "C. Module partition") +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")

save_panel(p4a, "Figure4A_corrected_full_data_partition_candidate", 4.8, 3.4)
save_panel(p4b, "Figure4B_corrected_full_data_coefficients_candidate", 5.2, 3.6)
save_panel(p4c, "Figure4C_corrected_full_data_modules_candidate", 7.4, 4.2)
save_combined(list(p4a, p4b, p4c), "Figure4_corrected_full_data_candidate")

cat("Regenerated combined and individual Figure 4 PDF/PNG/SVG files from Figure4 corrected full-data tables.\n")
