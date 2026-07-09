# Figure 3A endpoint-fix: signed gene-level screening plot.
# RStudio usage: open this file and click Source.
# Plotting only; no t-test is run.

PROJECT_DIR <- NULL
FDR_CUTOFF <- 0.01
SLOW_COLOR <- "#2166C2"
FAST_COLOR <- "#E8750B"
NS_COLOR <- "#A9A9A9"

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(path) && nzchar(path)) return(dirname(normalizePath(path)))
  }
  getwd()
}

base_dir <- if (is.null(PROJECT_DIR)) normalizePath(file.path(get_script_dir(), "..")) else PROJECT_DIR
input_path <- file.path(base_dir, "inputs", "endpointfix_Fig3A_screening_plot_table.tsv")
out_dir <- file.path(base_dir, "outputs")
panel_dir <- file.path(out_dir, "Figure3A")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
if (!requireNamespace("ggrepel", quietly = TRUE)) stop("Package 'ggrepel' is required.")
if (!requireNamespace("svglite", quietly = TRUE)) stop("Package 'svglite' is required.")

save_pdf_png_svg <- function(plot, file_stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(file_stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(file_stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(file_stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
}

message("Reading endpoint-fix Figure 3A input: ", input_path)
x <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("gene", "tvalue", "minus_log10_FDR", "FDR0.01_significant_old_rule", "direction_label")
missing_cols <- setdiff(required, names(x))
if (length(missing_cols) > 0) stop("Missing required endpoint-fix columns: ", paste(missing_cols, collapse = ", "))

x$tvalue <- as.numeric(x$tvalue)
x$minus_log10_FDR <- as.numeric(x$minus_log10_FDR)
x$significant <- as.logical(x$FDR0.01_significant_old_rule)
x$direction <- "non-significant"
x$direction[x$direction_label == "marine_slow_positive_t"] <- "marine-slow"
x$direction[x$direction_label == "marine_fast_negative_t"] <- "marine-fast"
x$direction <- factor(x$direction, levels = c("marine-slow", "marine-fast", "non-significant"))

label_genes <- c("SCN11A", "NSD1", "SLC44A5", "MYH13", "TGM3", "CAPN12", "TPX2")
x$label_gene <- ifelse(x$gene %in% label_genes & x$significant, x$gene, "")
counts <- table(x$direction)
legend_labels <- c(
  "marine-slow" = paste0("marine-slow (n = ", unname(counts["marine-slow"]), ")"),
  "marine-fast" = paste0("marine-fast (n = ", unname(counts["marine-fast"]), ")"),
  "non-significant" = paste0("non-significant (n = ", unname(counts["non-significant"]), ")")
)

plot_data <- x
plot_data$plot_order <- ifelse(plot_data$direction == "non-significant", 1, 2)
plot_data <- plot_data[order(plot_data$plot_order), ]
label_data <- plot_data[plot_data$label_gene != "", ]

write.table(x, file.path(panel_dir, "Figure3A_endpointfix_plot_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = tvalue, y = minus_log10_FDR, color = direction)) +
  ggplot2::geom_point(alpha = 0.75, size = 2.4, stroke = 0) +
  ggplot2::geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey55", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  ggrepel::geom_text_repel(data = label_data, ggplot2::aes(label = label_gene), color = "black", size = 3.2, min.segment.length = 0, box.padding = 0.35, point.padding = 0.2, max.overlaps = Inf, show.legend = FALSE) +
  ggplot2::scale_color_manual(values = c("marine-slow" = SLOW_COLOR, "marine-fast" = FAST_COLOR, "non-significant" = NS_COLOR), labels = legend_labels, drop = FALSE, name = NULL) +
  ggplot2::labs(x = "Signed t-value", y = expression(-log[10]~"(FDR)"), title = "Signed gene-level screening plot") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12), legend.position = c(0.20, 0.86), legend.background = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), axis.title = ggplot2::element_text(size = 11), axis.text = ggplot2::element_text(size = 9))

save_pdf_png_svg(p, file.path(out_dir, "Figure3A_endpointfix"), width = 6.0, height = 4.2)
message("Wrote endpoint-fix Figure 3A panel to: ", out_dir)
