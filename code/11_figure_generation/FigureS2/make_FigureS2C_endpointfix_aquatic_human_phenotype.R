# Supplementary Figure S2 endpointfix aquaticC: Monarch phenotype/measurement enrichment of aquatic-slow genes
#
# Current endpoint-fix status:
#   aquatic-slow has one accepted STRING/Monarch phenotype/measurement term
#   (EFO:0004503 Hematological measurement). Draw it as a one-term panel so
#   Figure S2 can retain the same A/B/C visual structure as the old figure.
#
# Expected input:
#   inputs/Figure3C_aquatic_plot_table.tsv
#
# Required columns:
#   #category, term description, observed gene count, strength, signal,
#   false discovery rate, selected_for_supp_figure, theme
#
# This script is designed to run by clicking Source in RStudio.

PROJECT_DIR <- NULL
INPUT_FILE <- "inputs/Figure3C_aquatic_plot_table.tsv"
OUT_DIR <- "outputs"
X_AXIS <- "observed gene count"

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

save_pdf_svg_png <- function(plot, file_stem, width, height) {
  ggplot2::ggsave(paste0(file_stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(file_stem, ".png"), plot, width = width, height = height, dpi = 300)
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
required <- c("#category", "term description", "observed gene count", "strength", "signal", "false discovery rate", "selected_for_supp_figure", "theme")
missing_cols <- setdiff(required, names(x))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

x <- x[x[["#category"]] == "Monarch", ]
x[[X_AXIS]] <- as.numeric(x[[X_AXIS]])
x[["observed gene count"]] <- as.numeric(x[["observed gene count"]])
x[["false discovery rate"]] <- as.numeric(x[["false discovery rate"]])
x$neg_log10_FDR <- -log10(pmax(x[["false discovery rate"]], .Machine$double.xmin))

plot_data <- x[x$selected_for_supp_figure %in% c(TRUE, "TRUE", "true", "1"), ]
if (nrow(plot_data) == 0) {
  plot_data <- x[x[["accepted_monarch_phenotype_measurement"]] %in% c(TRUE, "TRUE", "true", "1") |
                   x[["term ID"]] %in% c("EFO:0004503"), , drop = FALSE]
  if (nrow(plot_data) == 0) {
    stop("No accepted STRING/Monarch phenotype/measurement term was found for Figure S2C.")
  }
  plot_data$single_term_panel <- TRUE
  message("FigureS2C has one accepted term; drawing a one-term panel rather than stopping.")
}

plot_data <- plot_data[order(plot_data[[X_AXIS]], decreasing = FALSE), ]
plot_data$display_term <- factor(plot_data[["term description"]], levels = plot_data[["term description"]])

p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[X_AXIS]], y = display_term)) +
  ggplot2::geom_segment(ggplot2::aes(x = 0, xend = .data[[X_AXIS]], yend = display_term), color = "grey82", linewidth = 0.45) +
  ggplot2::geom_point(ggplot2::aes(size = .data[["observed gene count"]], color = neg_log10_FDR), alpha = 0.95) +
  ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
  ggplot2::scale_color_gradient(low = "#A7C7E7", high = "#08306B", name = expression(-log[10]~"(FDR)")) +
  ggplot2::scale_size_continuous(range = c(4.5, 9), name = "Gene count") +
  ggplot2::guides(
    color = ggplot2::guide_colorbar(title.position = "top", barheight = grid::unit(0.95, "in"), barwidth = grid::unit(0.22, "in")),
    size = ggplot2::guide_legend(title.position = "top", override.aes = list(color = "grey25"))
  ) +
  ggplot2::labs(x = "Gene count", y = NULL, title = "Phenotype/measurement enrichment of aquatic-slow genes") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 10.5),
    axis.text.y = ggplot2::element_text(size = 9),
    legend.position = "right",
    legend.title = ggplot2::element_text(size = 8.5),
    legend.text = ggplot2::element_text(size = 8.5),
    plot.margin = ggplot2::margin(10, 12, 8, 8)
  )

save_pdf_svg_png(p, file.path(out_dir, "FigureS2C_endpointfix_aquatic_human_phenotype"), width = 7.2, height = 2.85)
write.table(
  data.frame(
    status = if (nrow(plot_data) == 1) "GENERATED_SINGLE_TERM_PANEL" else "GENERATED",
    n_terms = nrow(plot_data),
    term_id = paste(plot_data[["term ID"]], collapse = ";"),
    term_label = paste(plot_data[["term description"]], collapse = ";"),
    FDR = paste(plot_data[["false discovery rate"]], collapse = ";"),
    gene_count = paste(plot_data[["observed gene count"]], collapse = ";"),
    notes = "Figure S2C was generated from accepted STRING/Monarch phenotype/measurement terms; endpoint-fix aquatic-slow contains one accepted EFO term.",
    stringsAsFactors = FALSE
  ),
  file.path(out_dir, "FigureS2C_endpointfix_generation_status.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
message("Wrote endpointfix Figure S2C outputs to: ", out_dir)
