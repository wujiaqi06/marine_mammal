#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggforce)
  library(svglite)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
work_dir <- file.path(root_dir, "09_figures", "New_Figures_endpointfix", "Figure4", "corrected_full_data_old_design")
input_dir <- file.path(work_dir, "inputs")
out_dir <- file.path(work_dir, "Figure4A_outputs")
dir.create(input_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

coef_path <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment", "Figure4_corrected_full_data", "Figure4_corrected_coefficient_architecture.tsv")
module_path <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment", "Figure4_corrected_full_data", "Figure4_corrected_module_partition.tsv")

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

coef <- read_tsv(coef_path)
module_part <- read_tsv(module_path)

coef$lasso_group <- ifelse(coef$predictor_class == "shared", "shared",
                           ifelse(coef$predictor_class == "marine_specific", "marine_specific", "aquatic_specific"))
coef$selected_in_marine <- coef$marine_nonzero
coef$selected_in_aquatic <- coef$aquatic_nonzero
coef$marine_coef <- coef$coefficient_marine
coef$aquatic_coef <- coef$coefficient_aquatic
coef$functional_annotation_short <- ifelse(coef$annotation_status == "REVIEW_UNANNOTATED", "not module-annotated", coef$candidate_module)
coef$functional_annotation_source <- coef$module_annotation_source
coef$annotation_confidence <- ifelse(coef$annotation_status == "REVIEW_UNANNOTATED", "unannotated", "module_annotated")
coef$notes <- ifelse(coef$annotation_status == "REVIEW_UNANNOTATED",
                     "Selected predictor retained; not assigned to a curated module for Figure 4C.",
                     "Selected predictor assigned to module annotation used for Figure 4C.")
coef$coefficient_direction_marine <- ifelse(coef$marine_coef > 0, "positive", ifelse(coef$marine_coef < 0, "negative", "zero"))
coef$coefficient_direction_aquatic <- ifelse(coef$aquatic_coef > 0, "positive", ifelse(coef$aquatic_coef < 0, "negative", "zero"))

rank_selected <- function(values, selected) {
  out <- rep(NA_integer_, length(values))
  idx <- which(selected)
  out[idx] <- rank(-abs(values[idx]), ties.method = "first")
  out
}
coef$marine_rank <- rank_selected(coef$marine_coef, coef$selected_in_marine)
coef$aquatic_rank <- rank_selected(coef$aquatic_coef, coef$selected_in_aquatic)

gene_level <- coef[, c(
  "gene", "lasso_group", "marine_coef", "aquatic_coef", "abs_marine_coef", "abs_aquatic_coef",
  "marine_rank", "aquatic_rank", "selected_in_marine", "selected_in_aquatic",
  "coefficient_direction_marine", "coefficient_direction_aquatic",
  "functional_annotation_short", "functional_annotation_source", "candidate_module",
  "annotation_confidence", "notes"
)]
write_tsv(gene_level, file.path(input_dir, "Table1_gene_level_working_table.corrected.tsv"))
write_tsv(gene_level, file.path(input_dir, "Table1_gene_level_working_table.tsv"))

marine_set <- sort(coef$gene[coef$selected_in_marine])
aquatic_set <- sort(coef$gene[coef$selected_in_aquatic])
shared <- intersect(marine_set, aquatic_set)
marine_only <- setdiff(marine_set, aquatic_set)
aquatic_only <- setdiff(aquatic_set, marine_set)

arch <- data.frame(
  comparison = "fix_marine_binary_vs_fix_aquatic_v2_corrected",
  marine_genes = length(marine_set),
  aquatic_genes = length(aquatic_set),
  shared_genes = length(shared),
  marine_only = length(marine_only),
  aquatic_only = length(aquatic_only),
  jaccard = length(shared) / length(union(marine_set, aquatic_set)),
  shared_marine_dominant = sum(coef$gene %in% shared & coef$abs_marine_coef > coef$abs_aquatic_coef),
  shared_aquatic_dominant = sum(coef$gene %in% shared & coef$abs_aquatic_coef > coef$abs_marine_coef),
  shared_balanced = sum(coef$gene %in% shared & coef$abs_aquatic_coef == coef$abs_marine_coef),
  stringsAsFactors = FALSE
)
write_tsv(arch, file.path(input_dir, "main_pair_architecture.tsv"))

display_to_col <- c(
  "Shared predictors" = "shared",
  "Marine-specific predictors" = "marine_specific",
  "Aquatic-specific predictors" = "aquatic_specific"
)
annotated <- module_part[module_part$candidate_module != "REVIEW_UNANNOTATED", , drop = FALSE]
modules <- unique(annotated$candidate_module)
module_summary <- data.frame(module = modules, stringsAsFactors = FALSE)
for (prefix in c("shared", "marine_specific", "aquatic_specific")) {
  module_summary[[paste0(prefix, "_genes")]] <- ""
  module_summary[[paste0("n_", prefix)]] <- 0L
  module_summary[[paste0("representative_", prefix, "_genes")]] <- ""
}
for (i in seq_len(nrow(annotated))) {
  category <- display_to_col[[annotated$display_category[[i]]]]
  row <- match(annotated$candidate_module[[i]], module_summary$module)
  module_summary[[paste0(category, "_genes")]][row] <- annotated$genes[[i]]
  module_summary[[paste0("n_", category)]][row] <- annotated$gene_count[[i]]
  genes <- strsplit(annotated$genes[[i]], ";", fixed = TRUE)[[1]]
  scores <- coef$max_abs_coef[match(genes, coef$gene)]
  genes <- genes[order(-scores, genes, na.last = TRUE)]
  module_summary[[paste0("representative_", category, "_genes")]][row] <- paste(head(genes, 5), collapse = "; ")
}

module_order <- c(
  "epithelial barrier / keratinization / body-surface interface",
  "hematological / platelet / circulatory regulation",
  "immune / cytokine / inflammatory regulation",
  "sensory remodeling / olfactory or visual systems",
  "ion channels / mechanosensation / membrane signaling",
  "sperm / cilia / CatSper / reproductive cell function",
  "DNA repair / chromatin / genome maintenance",
  "Skeletal, muscle and extracellular-matrix remodeling",
  "metabolism / mitochondrial / redox regulation",
  "fatty-acid / lipid metabolism",
  "endocrine / reproductive hormone systems",
  "transcriptional / systemic regulatory genes"
)
module_summary <- module_summary[match(module_order[module_order %in% module_summary$module], module_summary$module), , drop = FALSE]
module_summary$evidence_summary <- "Corrected Phase 12A full-data LASSO selected predictor module assignment; REVIEW_UNANNOTATED genes excluded from Figure 4C module rows."
module_summary$annotation_confidence <- "module_annotated"
module_summary$candidate_for_main_table <- "yes"
write_tsv(module_summary, file.path(input_dir, "Table1_module_summary_working.tsv"))
write_tsv(module_summary, file.path(input_dir, "Table1_main_text_candidate.tsv"))

unannotated <- module_part[module_part$candidate_module == "REVIEW_UNANNOTATED", c("display_category", "gene_count", "genes"), drop = FALSE]
write_tsv(unannotated, file.path(input_dir, "Figure4_unannotated_predictor_counts.tsv"))

plot_table <- data.frame(
  category = c("Marine-specific predictors", "Shared aquatic foundation", "Aquatic-specific predictors"),
  count = c(length(marine_only), length(shared), length(aquatic_only)),
  x = c(-0.88, 0, 0.88),
  y = c(0, 0, 0),
  label = c("Marine-specific\npredictors", "Shared aquatic\nfoundation", "Aquatic-specific\npredictors"),
  stringsAsFactors = FALSE
)
write_tsv(plot_table, file.path(out_dir, "Figure4A_plot_table.tsv"))

circle_table <- data.frame(
  set = c("Marine model", "Binary aquatic-dependence model"),
  x0 = c(-0.55, 0.55),
  y0 = c(0, 0),
  r = c(1.0, 1.0),
  stringsAsFactors = FALSE
)

p <- ggplot() +
  geom_circle(
    data = circle_table,
    aes(x0 = x0, y0 = y0, r = r, fill = set, color = set),
    alpha = 0.28,
    linewidth = 0.9
  ) +
  geom_text(data = plot_table, aes(x = x, y = y + 0.08, label = count),
            size = 7, fontface = "bold", color = "#111111") +
  geom_text(data = plot_table, aes(x = x, y = y - 0.22, label = label),
            size = 3.2, lineheight = 0.95, color = "#111111") +
  annotate("text", x = -0.78, y = 1.18, label = paste0("Marine model\n(n = ", length(marine_set), ")"),
           color = "#1F4E79", fontface = "bold", size = 3.5, lineheight = 0.95) +
  annotate("text", x = 0.78, y = 1.18, label = paste0("Binary aquatic-\ndependence model\n(n = ", length(aquatic_set), ")"),
           color = "#006D5B", fontface = "bold", size = 3.15, lineheight = 0.92) +
  annotate("text", x = 0, y = -1.24, label = paste0("Total baseline LASSO-selected predictors = ", length(union(marine_set, aquatic_set))),
           size = 3.3, color = "#111111") +
  scale_fill_manual(values = c("Marine model" = "#4C78A8", "Binary aquatic-dependence model" = "#44A08D")) +
  scale_color_manual(values = c("Marine model" = "#1F4E79", "Binary aquatic-dependence model" = "#006D5B")) +
  coord_equal(xlim = c(-1.85, 1.85), ylim = c(-1.45, 1.45), expand = FALSE) +
  labs(title = "LASSO gene-set partition") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0, size = 12, margin = margin(b = 8)),
    plot.margin = margin(8, 8, 8, 8)
  )

ggsave(file.path(out_dir, "Figure4A_gene_set_partition.pdf"), p, width = 4.7, height = 4.1, units = "in", useDingbats = FALSE, bg = "white")
ggsave(file.path(out_dir, "Figure4A_gene_set_partition.svg"), p, width = 4.7, height = 4.1, units = "in", bg = "white")
ggsave(file.path(out_dir, "Figure4A_gene_set_partition.png"), p, width = 4.7, height = 4.1, units = "in", dpi = 300, bg = "white")

cat("Corrected Figure 4A and inputs written to ", work_dir, "\n", sep = "")
