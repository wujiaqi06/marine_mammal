# Supplementary Figure S2 endpointfix aquaticA: signed gene-level screening plot
#
# Expected input:
#   inputs/Figure3A_aquatic_plot_table.tsv
#
# Required columns:
#   gene, tvalue, FDR_BH, neg_log10_FDR, direction, label_gene
#
# This script is designed to run by clicking Source in RStudio.

PROJECT_DIR <- NULL
INPUT_FILE <- "inputs/Figure3A_aquatic_plot_table.tsv"
OUT_DIR <- "outputs"
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

save_pdf_svg <- function(plot, file_stem, width, height) {
  ggplot2::ggsave(paste0(file_stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  if (requireNamespace("svglite", quietly = TRUE)) {
    ggplot2::ggsave(paste0(file_stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
  } else {
    warning("Package 'svglite' is not installed; SVG output was skipped for ", basename(file_stem), ".")
  }
}

base_dir <- if (is.null(PROJECT_DIR)) normalizePath(file.path(get_script_dir(), "..")) else PROJECT_DIR
input_path <- file.path(base_dir, INPUT_FILE)
out_dir <- file.path(base_dir, OUT_DIR)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

x <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("gene", "tvalue", "FDR_BH", "neg_log10_FDR", "direction", "label_gene")
missing_cols <- setdiff(required, names(x))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

x$direction <- factor(x$direction, levels = c("aquatic-slow", "aquatic-fast", "non-significant"))
counts <- table(x$direction)

plot_data <- x
plot_data$plot_order <- ifelse(plot_data$direction == "non-significant", 1, 2)
plot_data <- plot_data[order(plot_data$plot_order), ]
label_data <- plot_data[plot_data$label_gene != "", ]

legend_labels <- c(
  "aquatic-slow" = paste0("aquatic-slow (n = ", counts["aquatic-slow"], ")"),
  "aquatic-fast" = paste0("aquatic-fast (n = ", counts["aquatic-fast"], ")"),
  "non-significant" = paste0("non-significant (n = ", counts["non-significant"], ")")
)

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = tvalue, y = neg_log10_FDR, color = direction)) +
  ggplot2::geom_point(alpha = 0.75, size = 3, stroke = 0) +
  ggplot2::geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey55", linewidth = 0.35) +
  ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  ggplot2::scale_color_manual(
    values = c("aquatic-slow" = SLOW_COLOR, "aquatic-fast" = FAST_COLOR, "non-significant" = NS_COLOR),
    labels = legend_labels,
    drop = FALSE,
    name = NULL
  ) +
  ggplot2::labs(
    x = expression(italic(t)~"value"),
    y = expression(-log[10]~"(FDR)"),
    title = "Aquatic-dependence signed gene-level screening"
  ) +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    legend.position = c(0.18, 0.86),
    legend.background = ggplot2::element_blank(),
    legend.key = ggplot2::element_blank()
  )

if (nrow(label_data) > 0) {
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p + ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = label_gene),
      color = "black",
      size = 3.0,
      min.segment.length = 0,
      box.padding = 0.35,
      point.padding = 0.2,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  } else {
    p <- p + ggplot2::geom_text(data = label_data, ggplot2::aes(label = label_gene), color = "black", size = 3, hjust = -0.1, vjust = -0.3, show.legend = FALSE)
  }
}

save_pdf_svg(p, file.path(out_dir, "FigureS2A_endpointfix_aquatic_volcano"), width = 6.0, height = 4.2)
message("Wrote endpointfix Figure S2A outputs to: ", out_dir)

