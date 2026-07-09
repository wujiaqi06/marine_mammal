#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
})

root_dir <- normalizePath(Sys.getenv("MARINE_MAMMAL_ENDPOINTFIX_ROOT", unset = "."), mustWork = TRUE)
base_dir <- file.path(root_dir, "10_reviewer_risk_controls", "13_Figure4_Figure5_evidence_alignment")
fig_dir <- file.path(base_dir, "figures_draft")
fig5b_dir <- file.path(base_dir, "Figure5B_permutation_cleanup")
fig5c_dir <- file.path(base_dir, "Figure5C_corrected_turnover_module_null")
decision_dir <- file.path(base_dir, "decision")
qc_dir <- file.path(base_dir, "qc")
fig4_dir <- file.path(base_dir, "Figure4_corrected_full_data")

read_tsv <- function(path) read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")

save_all <- function(plot, stem, width, height) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height, units = "in", device = grDevices::pdf)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height, units = "in", dpi = 300)
  svglite::svglite(file.path(fig_dir, paste0(stem, ".svg")), width = width, height = height)
  print(plot)
  grDevices::dev.off()
}

theme_phase13 <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      legend.position = "bottom"
    )
}

## Patch 1: Result 2.3 neutral forward-link wording.
writeLines(c(
  "# Result 2.3 Forward-Link Patch",
  "",
  "Replace the old forward link about sparse predictor behaviour with:",
  "",
  "> This distributed single-gene asymmetry contrasts with the clade-dependent attenuation of sparse predictor models examined below.",
  "",
  "Do not use old wording that implies a strict null model unless explicitly describing deprecated legacy outputs. Under the Phase 13 nested sensitivity evidence layer, the marine drop-cetacean run is not strict-null."
), file.path(decision_dir, "Result2.3_forward_link_patch.md"))

## Patch 2: Figure 5B public-facing terminal-label counts only.
fig5b <- read_tsv(file.path(fig5b_dir, "Figure5B_endpointfix_source_check.tsv"))
fig5b$public_figure_count_scope <- "terminal_label_matched_permutation_design"
fig5b$branch_state_counts_public_display <- "no_internal_QC_only"
fig5b$branch_state_count_note <- "post-ASR branch-state counts are retained for QC only and should not appear in Figure 5B figure text, caption, or main Results."
write_tsv(fig5b, file.path(fig5b_dir, "Figure5B_endpointfix_source_check.tsv"))

writeLines(c(
  "# Figure 5B Caption Patch",
  "",
  "Figure 5B uses exact positive-count permutations for the endpoint-fix drop-cetacean marine single-gene screen.",
  "",
  "Public figure/caption text should use only the terminal-label matched permutation design:",
  "",
  "- Exact positive-count permutations",
  "- Matched terminal labels: 17 positives / 251 negatives",
  "- Observed: 894 FDR genes; 98.0% slow",
  "- n permutations = 200",
  "",
  "Observed endpoint-fix values: 894 FDR-significant genes, 876 slow and 18 fast (98.0% slow).",
  "",
  "Do not place post-ASR branch-state counts in the Figure 5B figure text, caption, or main Results. They remain only in `Figure5B_endpointfix_source_check.tsv` as internal QC."
), file.path(fig5b_dir, "Figure5B_caption_patch.md"))

p5b_text <- paste(
  "Exact positive-count permutations",
  "Matched terminal labels: 17 positives / 251 negatives",
  "Observed: 894 FDR genes; 98.0% slow",
  "n permutations = 200",
  sep = "\n"
)
p5b <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = "#F7F7F7", color = "#333333", linewidth = 0.4) +
  annotate("text", x = 0.05, y = 0.84, label = "Figure 5B endpoint-fix permutation control", hjust = 0, fontface = "bold", size = 4.2) +
  annotate("text", x = 0.05, y = 0.60, label = p5b_text, hjust = 0, lineheight = 1.15, size = 3.6) +
  annotate("text", x = 0.05, y = 0.12, label = "Branch-state counts are QC-only and not shown in public figure/caption text.", hjust = 0, size = 2.7, color = "#555555") +
  xlim(0, 1) + ylim(0, 1) +
  theme_void()
save_all(p5b, "Figure5B_endpointfix_cleaned_candidate", 5.2, 3.6)

## Patch 3: Figure 5C comparison-type audit and two candidate structures.
turnover <- read_tsv(file.path(fig5c_dir, "Figure5C_turnover_metrics.tsv"))
null_summary <- read_tsv(file.path(fig5c_dir, "Figure5C_module_null_summary.tsv"))

comparison_audit <- data.frame(
  comparison_id = c("5C_1", "5C_2", "5C_3"),
  set_A_name = c("marine baseline LASSO predictors", "marine baseline LASSO predictors", "binary aquatic baseline LASSO predictors"),
  set_B_name = c("whale-only LASSO predictors", "drop-cetacean marine-slow genes", "aquatic no-Cetacea LASSO predictors"),
  set_A_evidence_layer = c("Phase 12A corrected full-data LASSO predictors", "Phase 12A corrected full-data LASSO predictors", "Phase 12A corrected full-data LASSO predictors"),
  set_B_evidence_layer = c("Phase 12A corrected full-data LASSO predictors", "global branch-level t-test/FDR slow genes", "Phase 12A corrected full-data LASSO predictors"),
  predictor_vs_predictor = c("yes", "no", "yes"),
  cross_layer = c("no", "yes", "no"),
  recommended_main_or_supplement = c("main_candidate_2col_or_3col", "supplement_or_labeled_context", "main_candidate_2col_or_3col"),
  notes = c(
    "Same evidence layer; predictor-turnover comparison.",
    "Cross-layer comparison; descriptive context only and not equivalent to predictor-vs-predictor turnover.",
    "Same evidence layer; predictor-turnover comparison."
  ),
  stringsAsFactors = FALSE
)
write_tsv(comparison_audit, file.path(fig5c_dir, "Figure5C_comparison_type_audit.tsv"))

plot_fig5c_candidate <- function(include_cross_layer, stem, title, subtitle) {
  keep <- if (include_cross_layer) {
    turnover$comparison
  } else {
    turnover$comparison[turnover$comparison != "marine baseline vs drop-cetaceans slow genes"]
  }
  ns <- null_summary[null_summary$comparison %in% keep, ]
  ns$comparison_label <- ns$comparison
  ns$comparison_label[ns$comparison == "marine baseline vs whale-only"] <- "Marine baseline predictors\nvs whale-only predictors"
  ns$comparison_label[ns$comparison == "binary aquatic baseline vs aquatic no Cetacea"] <- "Aquatic baseline predictors\nvs aquatic no-Cetacea predictors"
  ns$comparison_label[ns$comparison == "marine baseline vs drop-cetaceans slow genes"] <- "Cross-layer context:\nmarine predictors vs drop-cetacean slow genes"
  ns$metric_label <- ns$metric
  ns$metric_label[ns$metric == "gene_overlap_Jaccard"] <- "Gene overlap\nJaccard"
  ns$metric_label[ns$metric == "module_presence_Jaccard"] <- "Module presence\nJaccard"
  ns$metric_label[ns$metric == "module_count_cosine_similarity"] <- "Module count\ncosine"
  ns$comparison_label <- factor(ns$comparison_label, levels = unique(ns$comparison_label))
  ns$metric_label <- factor(ns$metric_label, levels = c("Gene overlap\nJaccard", "Module presence\nJaccard", "Module count\ncosine"))
  p <- ggplot(ns, aes(metric_label, observed_value)) +
    geom_col(fill = "#5C6F82", width = 0.62, color = "grey25", linewidth = 0.2) +
    geom_point(aes(y = null_median), color = "#C74832", size = 2.1) +
    geom_text(aes(label = sprintf("%.2f", observed_value)), vjust = -0.35, size = 2.6) +
    facet_wrap(~ comparison_label, nrow = 1) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
      title = title,
      subtitle = subtitle,
      x = NULL,
      y = "Observed metric; red point = null median"
    ) +
    theme_phase13(8) +
    theme(strip.text = element_text(face = "bold"), legend.position = "none")
  if (include_cross_layer) {
    p <- p + labs(caption = "Middle column is cross-layer descriptive context, not the same class of predictor-turnover contrast.")
  }
  save_all(p, stem, ifelse(include_cross_layer, 9.2, 6.8), 4.0)
}

plot_fig5c_candidate(
  FALSE,
  "Figure5C_candidate_2col_predictor_vs_predictor",
  "Figure 5C candidate: predictor-vs-predictor contrasts only",
  "Cross-layer drop-cetacean slow-gene comparison moved to Supplementary/QC."
)
plot_fig5c_candidate(
  TRUE,
  "Figure5C_candidate_3col_with_cross_layer_context",
  "Figure 5C candidate: turnover with cross-layer context",
  "The middle column is explicitly labeled as a cross-layer descriptive comparison."
)

## Patch 4: Figure 4 caption/evidence note.
writeLines(c(
  "# Figure 4 Caption / Evidence Note Patch",
  "",
  "Figure 4 summarizes full-data final LASSO architecture after corrected preprocessing and is not a cross-validated validation panel.",
  "",
  "Module grouping is used as a descriptive organization of selected predictors and does not by itself imply module-level evolutionary recurrence; module recurrence across sensitivity comparisons is evaluated separately in Fig. 5C."
), file.path(fig4_dir, "Figure4_caption_evidence_note_patch.md"))

## Pro decision question for cross-layer column.
writeLines(c(
  "# Figure 5C Cross-Layer Column Decision Question for Pro",
  "",
  "Should Figure 5C main panel keep the cross-layer comparison \"marine baseline LASSO predictors vs drop-cetacean marine-slow genes,\" or should the main panel be restricted to predictor-vs-predictor comparisons and move the cross-layer comparison to Supplementary/QC?",
  "",
  "Control-chat concern: the cross-layer comparison has different evidence-layer semantics and may not be directly comparable to the predictor-vs-predictor columns. Its gene-overlap direction also differs from the two predictor-vs-predictor comparisons and may be dominated by set-size/layer differences."
), file.path(decision_dir, "Figure5C_cross_layer_column_for_Pro.md"))

## Update decision and run-summary notes without changing numerical outputs.
writeLines(c(
  "# Figure 5C Decision Memo",
  "",
  "Decision: Figure 5C requires Pro choice between two presentation structures.",
  "",
  "Candidate 5C-main-2col uses only predictor-vs-predictor comparisons:",
  "",
  "- Marine baseline predictors vs whale-only predictors.",
  "- Binary aquatic baseline predictors vs aquatic no-Cetacea predictors.",
  "",
  "Candidate 5C-main-3col retains the cross-layer comparison, but labels it explicitly as: cross-layer comparison, marine baseline LASSO predictors vs drop-cetacean marine-slow genes.",
  "",
  "The cross-layer comparison is descriptive context, not the same class of predictor-turnover contrast. Under the comparison-specific candidate-gene null, Figure 5C can support descriptive gene-level turnover, but should not be used for a strong above-null module recurrence claim."
), file.path(decision_dir, "Figure5C_turnover_module_null_decision_memo.md"))

writeLines(c(
  "# Figure 4-5 Alignment Decision Memo",
  "",
  "Overall classification: **ARCHITECTURE_REQUIRES_REVISION**.",
  "",
  "Figure 4 can be rebuilt from Phase 12A corrected full-data coefficients. The corrected full-data baseline architecture has 71 marine predictors, 101 binary aquatic-dependence predictors, 24 shared predictors, 47 marine-specific predictors, 77 aquatic-specific predictors, and 148 predictors in the union.",
  "",
  "Figure 5A can be rebuilt from nested supervised t-test feature-selection sensitivity values, but the previous strict-null drop-cetacean story changes. The nested `fix_drop_whale` run has AUC = 0.784626 and is classified as retained sparse prediction.",
  "",
  "Figure 5B is unaffected by LASSO nesting and was cleaned as endpoint-fix positive-count permutation provenance. Public figure/caption wording should use only terminal-label counts: 17 positives / 251 negatives, 894 FDR genes, and 98.0% slow.",
  "",
  "Figure 5C uses Figure 4-tied corrected full-data predictor sets and comparison-specific candidate-gene nulls. The package now provides both a 2-column predictor-vs-predictor candidate and a 3-column candidate that explicitly labels the cross-layer comparison.",
  "",
  "Implications:",
  "",
  "- Fig. 4: retain as corrected full-data architecture summary, with caption saying it is not a CV validation panel.",
  "- Fig. 5A: retain only with text rewritten around nested sensitivity values; do not use legacy or Phase 12A globalFDR AUCs as main.",
  "- Fig. 5B: retain with endpoint-fix terminal-label matched permutation wording.",
  "- Fig. 5C: ask Pro whether the cross-layer comparison belongs in main or should move to Supplementary/QC."
), file.path(decision_dir, "Figure4_Figure5_alignment_decision_memo.md"))

cat("Opus4.7 Phase 13 patch applied.\n")
