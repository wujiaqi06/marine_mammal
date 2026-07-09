# Figure 3B endpoint-fix: directional asymmetry across marine analyses.
# RStudio usage: open this file and click Source.
# Plotting only; no t-test is run.

PROJECT_DIR <- NULL
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

base_dir <- if (is.null(PROJECT_DIR)) normalizePath(file.path(get_script_dir(), "..")) else PROJECT_DIR
input_path <- file.path(base_dir, "inputs", "endpointfix_Fig3B_directional_asymmetry_summary.tsv")
out_dir <- file.path(base_dir, "outputs")
panel_dir <- file.path(out_dir, "Figure3B")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
if (!requireNamespace("svglite", quietly = TRUE)) stop("Package 'svglite' is required.")

save_pdf_png_svg <- function(plot, file_stem, width, height, dpi = 300) {
  ggplot2::ggsave(paste0(file_stem, ".pdf"), plot, width = width, height = height, device = "pdf", useDingbats = FALSE)
  ggplot2::ggsave(paste0(file_stem, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
  ggplot2::ggsave(paste0(file_stem, ".svg"), plot, width = width, height = height, device = svglite::svglite)
}

message("Reading endpoint-fix Figure 3B input: ", input_path)
x <- read.delim(input_path, check.names = FALSE, stringsAsFactors = FALSE)

run_map <- data.frame(
  run_id = c("fix_marine_binary", "fix_drop_whale", "fix_drop_polar_bear", "fix_drop_sea_otter", "fix_pinniped_only"),
  run_label = c("baseline", "drop\ncetaceans", "drop\npolar bear", "drop\nsea otter", "pinniped-\nonly"),
  run_label_plain = c("baseline", "drop cetaceans", "drop polar bear", "drop sea otter", "pinniped-only"),
  stringsAsFactors = FALSE
)
x <- merge(run_map, x, by = "run_id", all.x = TRUE, sort = FALSE)
x <- x[match(run_map$run_id, x$run_id), ]
if (anyNA(x$marine_slow_positive_t)) stop("Missing expected Figure 3B run in endpoint-fix table.")

x$marine_slow_positive_t <- as.numeric(x$marine_slow_positive_t)
x$marine_fast_negative_t <- as.numeric(x$marine_fast_negative_t)
x$slow_pct <- 100 * as.numeric(x$slow_proportion)
x$fast_pct <- 100 - x$slow_pct
x$slow_pct_label <- sprintf("%.1f%%", x$slow_pct)
x$run_label <- factor(x$run_label, levels = run_map$run_label)

plot_long <- rbind(
  data.frame(run_label = x$run_label, direction = "marine-fast", proportion = x$fast_pct),
  data.frame(run_label = x$run_label, direction = "marine-slow", proportion = x$slow_pct)
)
plot_long$direction <- factor(plot_long$direction, levels = c("marine-fast", "marine-slow"))

write.table(x, file.path(panel_dir, "Figure3B_endpointfix_plot_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(run_map, file.path(panel_dir, "Figure3B_endpointfix_run_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = run_label, y = proportion, fill = direction)) +
  ggplot2::geom_col(width = 0.62, color = "black", linewidth = 0.25) +
  ggplot2::geom_text(data = x, ggplot2::aes(x = run_label, y = 104, label = slow_pct_label), inherit.aes = FALSE, size = 3.7, fontface = "bold") +
  ggplot2::scale_fill_manual(values = c("marine-slow" = SLOW_COLOR, "marine-fast" = FAST_COLOR), breaks = c("marine-slow", "marine-fast"), name = NULL) +
  ggplot2::scale_y_continuous(limits = c(0, 110), breaks = c(0, 25, 50, 75, 100), expand = ggplot2::expansion(mult = c(0, 0.02))) +
  ggplot2::labs(x = NULL, y = "Proportion of significant genes (%)", title = "Directional asymmetry across analyses") +
  ggplot2::theme_classic(base_size = 11) +
  ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 12), axis.text.x = ggplot2::element_text(size = 9), axis.title.y = ggplot2::element_text(size = 11), legend.position = "bottom", legend.key.size = grid::unit(0.45, "cm"))

save_pdf_png_svg(p, file.path(out_dir, "Figure3B_endpointfix"), width = 5.0, height = 4.0)
message("Wrote endpoint-fix Figure 3B panel to: ", out_dir)
