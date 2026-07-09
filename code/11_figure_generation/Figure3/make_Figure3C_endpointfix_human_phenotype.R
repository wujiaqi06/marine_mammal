# Figure 3C endpoint-fix: phenotype and measurement enrichment of marine-slow genes.
#
# RStudio usage:
#   Open this file and click Source.
#
# This endpoint-fix script reads:
#   ../inputs/endpointfix_Fig3C_enrichment_terms_for_plot.tsv
#
# It does not run enrichment. It only plots the selected endpoint-fix terms.

PROJECT_DIR <- NULL

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
  }
  if (requireNamespace("rstudioapi", quietly = TRUE)) {
    path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(path) && nzchar(path)) {
      return(dirname(normalizePath(path)))
    }
  }
  getwd()
}

base_dir <- if (is.null(PROJECT_DIR)) normalizePath(file.path(get_script_dir(), "..")) else PROJECT_DIR
input_path <- file.path(base_dir, "inputs", "endpointfix_Fig3C_enrichment_terms_for_plot.tsv")
out_dir <- file.path(base_dir, "outputs")
panel_dir <- file.path(out_dir, "Figure3C")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
if (!requireNamespace("svglite", quietly = TRUE)) stop("Package 'svglite' is required.")

save_pdf_png_svg <- function(plot, file_stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(file_stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(file_stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(file_stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
}

message("Reading endpoint-fix Figure 3C input: ", input_path)
x <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)

required <- c("gene_set", "category", "term_id", "term_description", "observed_gene_count",
              "FDR", "minus_log10_FDR", "selected_for_main_figure")
missing_cols <- setdiff(required, names(x))
if (length(missing_cols) > 0) {
  stop("Missing required endpoint-fix columns: ", paste(missing_cols, collapse = ", "))
}

x$observed_gene_count <- as.numeric(x$observed_gene_count)
x$FDR <- as.numeric(x$FDR)
x$minus_log10_FDR <- as.numeric(x$minus_log10_FDR)
x$selected_for_main_figure <- as.logical(x$selected_for_main_figure)

plot_data <- x[x$gene_set == "marine_slow" & x$selected_for_main_figure, ]
if (nrow(plot_data) == 0) {
  stop("No selected endpoint-fix terms available for Figure 3C.")
}

plot_data <- plot_data[order(plot_data$observed_gene_count, decreasing = FALSE), ]
plot_data$display_term <- factor(plot_data$term_description, levels = plot_data$term_description)

write.table(
  x,
  file.path(panel_dir, "Figure3C_endpointfix_full_enrichment_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
write.table(
  plot_data,
  file.path(panel_dir, "Figure3C_endpointfix_plot_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = observed_gene_count, y = display_term)) +
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

save_pdf_png_svg(p, file.path(out_dir, "Figure3C_endpointfix"), width = 7.2, height = 4.5)
message("Wrote endpoint-fix Figure 3C panel to: ", out_dir)
