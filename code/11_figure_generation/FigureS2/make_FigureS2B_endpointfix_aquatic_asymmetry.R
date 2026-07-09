# Supplementary Figure S2 endpointfix aquaticB: directional asymmetry across aquatic_v2 drop analyses
#
# Expected input:
#   inputs/Figure3B_aquatic_plot_table.tsv
#
# Required columns:
#   run_id, run_label, aquatic_slow_n, aquatic_fast_n, aquatic_slow_pct
#
# This script is designed to run by clicking Source in RStudio.

PROJECT_DIR <- NULL
INPUT_FILE <- "inputs/Figure3B_aquatic_plot_table.tsv"
OUT_DIR <- "outputs"

SLOW_COLOR <- "#2166C2"
FAST_COLOR <- "#E8750B"

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
required <- c("run_id", "run_label", "aquatic_slow_n", "aquatic_fast_n", "aquatic_slow_pct")
missing_cols <- setdiff(required, names(x))
if (length(missing_cols) > 0) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

run_order <- x$run_label
x$run_label <- factor(x$run_label, levels = run_order)
x$aquatic_slow_n <- as.numeric(x$aquatic_slow_n)
x$aquatic_fast_n <- as.numeric(x$aquatic_fast_n)
x$total_sig_n <- x$aquatic_slow_n + x$aquatic_fast_n
x$aquatic_slow_pct <- 100 * x$aquatic_slow_n / x$total_sig_n
x$aquatic_fast_pct <- 100 * x$aquatic_fast_n / x$total_sig_n
x$slow_pct_label <- sprintf("%.1f%%", x$aquatic_slow_pct)

plot_long <- rbind(
  data.frame(run_label = x$run_label, direction = "aquatic-fast", count = x$aquatic_fast_n, proportion = x$aquatic_fast_pct),
  data.frame(run_label = x$run_label, direction = "aquatic-slow", count = x$aquatic_slow_n, proportion = x$aquatic_slow_pct)
)
plot_long$direction <- factor(plot_long$direction, levels = c("aquatic-fast", "aquatic-slow"))

p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = run_label, y = proportion, fill = direction)) +
  ggplot2::geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  ggplot2::geom_text(data = x, ggplot2::aes(x = run_label, y = 103, label = slow_pct_label), inherit.aes = FALSE, size = 3.2, fontface = "bold") +
  ggplot2::scale_fill_manual(values = c("aquatic-slow" = SLOW_COLOR, "aquatic-fast" = FAST_COLOR), breaks = c("aquatic-slow", "aquatic-fast"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, 108), breaks = c(0, 25, 50, 75, 100), expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::labs(x = NULL, y = "Proportion of significant genes (%)", title = "Aquatic directional asymmetry across analyses") +
  ggplot2::theme_classic(base_size = 10.5) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold", size = 12),
    axis.text.x = ggplot2::element_text(size = 8, angle = 35, hjust = 1),
    legend.position = "bottom"
  )

save_pdf_svg(p, file.path(out_dir, "FigureS2B_endpointfix_aquatic_asymmetry"), width = 6.8, height = 4.2)
message("Wrote endpointfix Figure S2B outputs to: ", out_dir)
