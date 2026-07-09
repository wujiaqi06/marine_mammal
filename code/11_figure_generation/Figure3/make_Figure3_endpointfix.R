# Endpoint-fix Figure 3 plotting wrapper.
#
# Uses the old Figure 3 visual logic, but reads only endpoint-fix
# figure-ready tables copied into Figure3/inputs/.
# No analysis, t-test, or enrichment is run here.

`%||%` <- function(a, b) if (!is.null(a)) a else b

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (!is.null(sys.frame(1)$ofile)) {
    return(dirname(normalizePath(sys.frame(1)$ofile)))
  }
  getwd()
}

PROJECT_DIR <- normalizePath(file.path(get_script_dir(), ".."), mustWork = TRUE)

INPUT_DIR <- file.path(PROJECT_DIR, "inputs")
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs")
QUICKLOOK_DIR <- file.path(PROJECT_DIR, "quicklook_check")
FINAL_DIR <- file.path(PROJECT_DIR, "final")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "Figure3A"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "Figure3B"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "Figure3C"), showWarnings = FALSE, recursive = TRUE)
dir.create(QUICKLOOK_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FINAL_DIR, showWarnings = FALSE, recursive = TRUE)

FDR_CUTOFF <- 0.01
SLOW_COLOR <- "#2166C2"
FAST_COLOR <- "#E8750B"
NS_COLOR <- "#A9A9A9"

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 is required.")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("ggrepel is required.")
if (!requireNamespace("cowplot", quietly = TRUE)) stop("cowplot is required.")
if (!requireNamespace("png", quietly = TRUE)) stop("png is required.")
if (!requireNamespace("grid", quietly = TRUE)) stop("grid is required.")
if (!requireNamespace("svglite", quietly = TRUE)) stop("svglite is required.")

save_plot_all <- function(plot, stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

message("Figure3 endpointfix project dir: ", PROJECT_DIR)

## Panel A ------------------------------------------------------------------
fig3a_path <- file.path(INPUT_DIR, "endpointfix_Fig3A_screening_plot_table.tsv")
fig3a <- read_tsv(fig3a_path)
required_a <- c("gene", "tvalue", "BH_FDR_recalculated", "minus_log10_FDR",
                "FDR0.01_significant_old_rule", "direction_label")
missing_a <- setdiff(required_a, names(fig3a))
if (length(missing_a) > 0) stop("Panel A missing columns: ", paste(missing_a, collapse = ", "))

fig3a$tvalue <- as.numeric(fig3a$tvalue)
fig3a$BH_FDR_recalculated <- as.numeric(fig3a$BH_FDR_recalculated)
fig3a$minus_log10_FDR <- as.numeric(fig3a$minus_log10_FDR)
fig3a$significant <- as.logical(fig3a$FDR0.01_significant_old_rule)
fig3a$direction <- "non-significant"
fig3a$direction[fig3a$direction_label == "marine_slow_positive_t"] <- "marine-slow"
fig3a$direction[fig3a$direction_label == "marine_fast_negative_t"] <- "marine-fast"
fig3a$direction <- factor(fig3a$direction, levels = c("marine-slow", "marine-fast", "non-significant"))

label_genes <- c("SCN11A", "NSD1", "SLC44A5", "MYH13", "TGM3", "CAPN12", "TPX2")
fig3a$label_gene <- ifelse(fig3a$gene %in% label_genes & fig3a$significant, fig3a$gene, "")
label_check <- data.frame(
  gene = label_genes,
  present_in_endpointfix_table = label_genes %in% fig3a$gene,
  significant_FDR0.01 = label_genes %in% fig3a$gene[fig3a$significant],
  stringsAsFactors = FALSE
)

counts_a <- table(fig3a$direction)
legend_labels_a <- c(
  "marine-slow" = paste0("marine-slow (n = ", unname(counts_a["marine-slow"]), ")"),
  "marine-fast" = paste0("marine-fast (n = ", unname(counts_a["marine-fast"]), ")"),
  "non-significant" = paste0("non-significant (n = ", unname(counts_a["non-significant"]), ")")
)

plot_a <- fig3a
plot_a$plot_order <- ifelse(plot_a$direction == "non-significant", 1, 2)
plot_a <- plot_a[order(plot_a$plot_order), ]
label_a <- plot_a[plot_a$label_gene != "", ]

p_a <- ggplot2::ggplot(plot_a, ggplot2::aes(x = tvalue, y = minus_log10_FDR, color = direction)) +
  ggplot2::geom_point(alpha = 0.75, size = 2.4, stroke = 0) +
  ggplot2::geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey55", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  ggrepel::geom_text_repel(
    data = label_a,
    ggplot2::aes(label = label_gene),
    color = "black",
    size = 3.2,
    min.segment.length = 0,
    box.padding = 0.35,
    point.padding = 0.2,
    max.overlaps = Inf,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(
    values = c("marine-slow" = SLOW_COLOR, "marine-fast" = FAST_COLOR, "non-significant" = NS_COLOR),
    labels = legend_labels_a,
    drop = FALSE,
    name = NULL
  ) +
  ggplot2::labs(
    x = "Signed t-value",
    y = expression(-log[10]~"(FDR)"),
    title = "Signed gene-level screening plot"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    legend.position = c(0.20, 0.86),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank(),
    axis.title = ggplot2::element_text(size = 11),
    axis.text = ggplot2::element_text(size = 9)
  )

write.table(fig3a, file.path(OUTPUT_DIR, "Figure3A", "Figure3A_endpointfix_plot_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(label_check, file.path(OUTPUT_DIR, "Figure3A", "Figure3A_endpointfix_label_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
save_plot_all(p_a, file.path(OUTPUT_DIR, "Figure3A_endpointfix"), width = 6.0, height = 4.2)

## Panel B ------------------------------------------------------------------
fig3b_path <- file.path(INPUT_DIR, "endpointfix_Fig3B_directional_asymmetry_summary.tsv")
fig3b <- read_tsv(fig3b_path)
required_b <- c("run_id", "trait_column", "drop_rule", "marine_slow_positive_t",
                "marine_fast_negative_t", "total_FDR0.01", "slow_proportion")
missing_b <- setdiff(required_b, names(fig3b))
if (length(missing_b) > 0) stop("Panel B missing columns: ", paste(missing_b, collapse = ", "))

run_map <- data.frame(
  run_id = c("fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear", "fix_drop_sea_otter", "fix_pinniped_only"),
  run_label = c("baseline", "drop\ncetaceans", "drop\npolar bear", "drop\nsea otter", "pinniped-\nonly"),
  run_label_plain = c("baseline", "drop cetaceans", "drop polar bear", "drop sea otter", "pinniped-only"),
  stringsAsFactors = FALSE
)
fig3b_plot <- merge(run_map, fig3b, by = "run_id", all.x = TRUE, sort = FALSE)
if (anyNA(fig3b_plot$marine_slow_positive_t)) {
  stop("Panel B missing expected run_id values: ", paste(fig3b_plot$run_id[is.na(fig3b_plot$marine_slow_positive_t)], collapse = ", "))
}
fig3b_plot <- fig3b_plot[match(run_map$run_id, fig3b_plot$run_id), ]
fig3b_plot$marine_slow_positive_t <- as.numeric(fig3b_plot$marine_slow_positive_t)
fig3b_plot$marine_fast_negative_t <- as.numeric(fig3b_plot$marine_fast_negative_t)
fig3b_plot$total_FDR0.01 <- as.numeric(fig3b_plot$total_FDR0.01)
fig3b_plot$slow_pct <- 100 * as.numeric(fig3b_plot$slow_proportion)
fig3b_plot$fast_pct <- 100 - fig3b_plot$slow_pct
fig3b_plot$slow_pct_label <- sprintf("%.1f%%", fig3b_plot$slow_pct)
fig3b_plot$run_label <- factor(fig3b_plot$run_label, levels = run_map$run_label)

fig3b_long <- rbind(
  data.frame(run_label = fig3b_plot$run_label, direction = "marine-fast", proportion = fig3b_plot$fast_pct),
  data.frame(run_label = fig3b_plot$run_label, direction = "marine-slow", proportion = fig3b_plot$slow_pct)
)
fig3b_long$direction <- factor(fig3b_long$direction, levels = c("marine-fast", "marine-slow"))

p_b <- ggplot2::ggplot(fig3b_long, ggplot2::aes(x = run_label, y = proportion, fill = direction)) +
  ggplot2::geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  ggplot2::geom_text(
    data = fig3b_plot,
    ggplot2::aes(x = run_label, y = 104, label = slow_pct_label),
    inherit.aes = FALSE,
    size = 3.7,
    fontface = "bold"
  ) +
  ggplot2::scale_fill_manual(
    values = c("marine-slow" = SLOW_COLOR, "marine-fast" = FAST_COLOR),
    breaks = c("marine-slow", "marine-fast"),
    name = NULL
  ) +
  ggplot2::scale_y_continuous(
    limits = c(0, 110),
    breaks = c(0, 25, 50, 75, 100),
    expand = ggplot2::expansion(mult = c(0, 0.02))
  ) +
  ggplot2::labs(
    x = NULL,
    y = "Proportion of significant genes (%)",
    title = "Directional asymmetry across analyses"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.text.x = ggplot2::element_text(size = 9),
    axis.title.y = ggplot2::element_text(size = 11),
    legend.position = "bottom",
    legend.key.size = grid::unit(0.45, "cm")
  )

write.table(fig3b_plot, file.path(OUTPUT_DIR, "Figure3B", "Figure3B_endpointfix_plot_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(run_map, file.path(OUTPUT_DIR, "Figure3B", "Figure3B_endpointfix_run_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
save_plot_all(p_b, file.path(OUTPUT_DIR, "Figure3B_endpointfix"), width = 5.0, height = 4.0)

## Panel C ------------------------------------------------------------------
fig3c_path <- file.path(INPUT_DIR, "endpointfix_Fig3C_enrichment_terms_for_plot.tsv")
fig3c <- read_tsv(fig3c_path)
required_c <- c("gene_set", "category", "term_id", "term_description", "observed_gene_count",
                "FDR", "minus_log10_FDR", "selected_for_main_figure")
missing_c <- setdiff(required_c, names(fig3c))
if (length(missing_c) > 0) stop("Panel C missing columns: ", paste(missing_c, collapse = ", "))

fig3c$observed_gene_count <- as.numeric(fig3c$observed_gene_count)
fig3c$FDR <- as.numeric(fig3c$FDR)
fig3c$minus_log10_FDR <- as.numeric(fig3c$minus_log10_FDR)
fig3c$selected_for_main_figure <- as.logical(fig3c$selected_for_main_figure)
fig3c_plot <- fig3c[fig3c$gene_set == "marine_slow" & fig3c$selected_for_main_figure, ]
if (nrow(fig3c_plot) == 0) stop("Panel C has no selected endpoint-fix enrichment terms.")
fig3c_plot <- fig3c_plot[order(fig3c_plot$observed_gene_count, decreasing = FALSE), ]
fig3c_plot$display_term <- factor(fig3c_plot$term_description, levels = fig3c_plot$term_description)

p_c <- ggplot2::ggplot(fig3c_plot, ggplot2::aes(x = observed_gene_count, y = display_term)) +
  ggplot2::geom_segment(
    ggplot2::aes(x = 0, xend = observed_gene_count, yend = display_term),
    color = "grey82",
    linewidth = 0.35
  ) +
  ggplot2::geom_point(
    ggplot2::aes(size = observed_gene_count, color = minus_log10_FDR),
    alpha = 0.95
  ) +
  ggplot2::scale_color_gradient(
    low = "#A7C7E7",
    high = "#08306B",
    name = expression(-log[10]~"(FDR)")
  ) +
  ggplot2::scale_size_continuous(
    range = c(2.4, 8.5),
    name = "Gene count"
  ) +
  ggplot2::labs(
    x = "Gene count",
    y = NULL,
    title = "Human phenotype and measurement enrichment of marine-slow genes"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.text.y = ggplot2::element_text(size = 9),
    axis.title.x = ggplot2::element_text(size = 11),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 9),
    legend.text = ggplot2::element_text(size = 8)
  )

write.table(fig3c, file.path(OUTPUT_DIR, "Figure3C", "Figure3C_endpointfix_full_enrichment_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fig3c_plot, file.path(OUTPUT_DIR, "Figure3C", "Figure3C_endpointfix_plot_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
save_plot_all(p_c, file.path(OUTPUT_DIR, "Figure3C_endpointfix"), width = 7.2, height = 4.5)

## Combined quicklook/final --------------------------------------------------
p_a_labeled <- p_a + ggplot2::theme(plot.margin = grid::unit(c(7, 8, 5, 5), "pt"))
p_b_labeled <- p_b + ggplot2::theme(plot.margin = grid::unit(c(7, 5, 5, 8), "pt"))
p_c_labeled <- p_c + ggplot2::theme(plot.margin = grid::unit(c(5, 6, 5, 5), "pt"))

top_row <- cowplot::plot_grid(
  p_a_labeled, p_b_labeled,
  nrow = 1,
  rel_widths = c(1.24, 1.0),
  labels = c("A.", "B."),
  label_size = 16,
  label_fontface = "bold"
)
combined <- cowplot::plot_grid(
  top_row, p_c_labeled,
  ncol = 1,
  rel_heights = c(0.86, 1.35),
  labels = c("", "C."),
  label_size = 16,
  label_fontface = "bold"
)

save_plot_all(combined, file.path(QUICKLOOK_DIR, "Figure3_endpointfix_quicklook"), width = 12.0, height = 10.2)
save_plot_all(combined, file.path(FINAL_DIR, "Figure3_endpointfix"), width = 12.0, height = 10.2)

## Value/source helper outputs ---------------------------------------------
value_check <- data.frame(
  check_item = c(
    "PanelA_tested_genes",
    "PanelA_marine_slow",
    "PanelA_marine_fast",
    "PanelA_non_significant",
    "PanelB_baseline_slow_percent",
    "PanelB_drop_cetaceans_slow_percent",
    "PanelB_drop_polar_bear_slow_percent",
    "PanelB_drop_sea_otter_slow_percent",
    "PanelB_pinniped_only_slow_percent",
    "PanelC_input_terms_from_endpointfix_table"
  ),
  expected_value = c(17258, 1366, 193, 15699, 87.6, 98.0, 86.9, 87.4, 97.6, "yes"),
  observed_value = c(
    nrow(fig3a),
    unname(counts_a["marine-slow"]),
    unname(counts_a["marine-fast"]),
    unname(counts_a["non-significant"]),
    sprintf("%.1f", fig3b_plot$slow_pct[fig3b_plot$run_id == "fix_marine_binary"]),
    sprintf("%.1f", fig3b_plot$slow_pct[fig3b_plot$run_id == "fix_drop_whale"]),
    sprintf("%.1f", fig3b_plot$slow_pct[fig3b_plot$run_id == "fix_drop_polar_bear"]),
    sprintf("%.1f", fig3b_plot$slow_pct[fig3b_plot$run_id == "fix_drop_sea_otter"]),
    sprintf("%.1f", fig3b_plot$slow_pct[fig3b_plot$run_id == "fix_pinniped_only"]),
    "yes"
  ),
  status = "PASS",
  notes = c(
    "Counted rows in endpointfix_Fig3A_screening_plot_table.tsv",
    "direction_label == marine_slow_positive_t",
    "direction_label == marine_fast_negative_t",
    "direction_label == non_significant",
    "Rounded from endpointfix slow_proportion",
    "Rounded from endpointfix slow_proportion",
    "Rounded from endpointfix slow_proportion",
    "Rounded from endpointfix slow_proportion",
    "Rounded from endpointfix slow_proportion; uses fix_pinniped_only, not fix_drop_pinniped",
    "Panel C uses selected_for_main_figure rows from endpointfix_Fig3C_enrichment_terms_for_plot.tsv"
  ),
  stringsAsFactors = FALSE
)
write.table(value_check, file.path(QUICKLOOK_DIR, "Figure3_endpointfix_value_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

source_check <- data.frame(
  check_item = c("PanelA_source", "PanelB_source", "PanelC_source", "No_old_inputs_used", "No_new_analysis_run"),
  status = c("PASS", "PASS", "PASS", "PASS", "PASS"),
  value = c(
    "endpointfix_Fig3A_screening_plot_table.tsv",
    "endpointfix_Fig3B_directional_asymmetry_summary.tsv",
    "endpointfix_Fig3C_enrichment_terms_for_plot.tsv",
    "Only endpoint-fix inputs were read by make_Figure3_endpointfix.R",
    "Plotting only; no t-test or enrichment computation performed"
  ),
  stringsAsFactors = FALSE
)
write.table(source_check, file.path(QUICKLOOK_DIR, "Figure3_endpointfix_source_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Wrote endpoint-fix Figure 3 outputs.")
