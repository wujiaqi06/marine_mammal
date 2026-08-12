#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(svglite)
})

args_full <- commandArgs(FALSE)
file_arg <- grep("^--file=", args_full, value = TRUE)
script_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "."
package_root <- normalizePath(file.path(dirname(script_file), "..", ".."), mustWork = TRUE)
source_root <- Sys.getenv(
  "MARINE_MAMMAL_FIGS3_SOURCE_ROOT",
  unset = file.path(package_root, "source_data", "FigS3_ancestor_fingerprints")
)
out_dir <- Sys.getenv(
  "MARINE_MAMMAL_FIGS3_OUTPUT",
  unset = file.path(package_root, "reproduction_outputs", "FigureS3")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fp_file <- file.path(source_root, "SourceData_FigS3_ancestor_fingerprints_long.csv")
profile_file <- file.path(source_root, "SourceData_FigS3_ancestor_profile_scores.csv")
if (!file.exists(fp_file) || !file.exists(profile_file)) {
  stop("Fig. S3 source tables are missing. Set MARINE_MAMMAL_FIGS3_SOURCE_ROOT if they are stored elsewhere.")
}

fp <- read.csv(fp_file, stringsAsFactors = FALSE, check.names = FALSE)
profiles <- read.csv(profile_file, stringsAsFactors = FALSE, check.names = FALSE)

INK <- "#1F2430"
MUTED <- "#5F6678"
GRID <- "#E6E8F0"
CARD <- "#D7DCE5"
RED <- "#8B0000"
BLUE_NEG <- "#1F4E79"

row_order <- c(
  "Cetacea + hippo ancestor",
  "Crown Cetacea",
  "Mysticeti",
  "Odontoceti",
  "Crown Pinnipedia",
  "Phocidae",
  "Otarioidea",
  "Otariidae",
  "Otariinae",
  "Sirenia",
  "Dugong + Hydrodamalis",
  "Brown bear + polar bear",
  "Lutrinae / Enhydra-related"
)
row_order <- row_order[row_order %in% unique(fp$display_label)]
fp <- fp[fp$display_label %in% row_order, , drop = FALSE]

profile_lookup <- profiles[, c("display_label", "branch_group", "marine_profile_score", "aquatic_profile_score"), drop = FALSE]
profile_lookup$branch_group_short <- ifelse(
  profile_lookup$branch_group == "Cetacean axis", "cetacean internal branch",
  ifelse(profile_lookup$branch_group == "Pinniped axis", "pinniped internal branch", "sirenian / edge internal branch")
)
fp <- merge(fp, profile_lookup, by = "display_label", all.x = TRUE, sort = FALSE)
fp$display_label <- factor(fp$display_label, levels = row_order)
fp$bar_sign <- ifelse(fp$contribution >= 0, "positive contribution", "negative contribution")

label_top_genes <- function(d, n = 6) {
  do.call(rbind, lapply(split(d, paste(d$display_label, d$model, sep = "||")), function(part) {
    part <- part[order(part$abs_contribution, decreasing = TRUE), , drop = FALSE]
    part$label_rank <- seq_len(nrow(part))
    part$show_label <- part$label_rank <= n
    part
  }))
}
fp <- label_top_genes(fp, n = 7)

ylim <- max(abs(fp$contribution), na.rm = TRUE) * 1.27
header_df <- unique(fp[, c("display_label", "branch_group_short", "marine_profile_score", "aquatic_profile_score", "model"), drop = FALSE])
header_df$profile_score <- ifelse(
  header_df$model == "marine",
  header_df$marine_profile_score,
  header_df$aquatic_profile_score
)
header_df$header_label <- paste0(
  as.character(header_df$display_label),
  ", profile p=", sprintf("%.3f", header_df$profile_score),
  "\n", header_df$branch_group_short
)
header_df$gene_order <- 1
header_df$contribution <- ylim * 0.80

plot_one_axis <- function(model_name, title, subtitle, n_pred) {
  d <- fp[fp$model == model_name, , drop = FALSE]
  h <- header_df[header_df$model == model_name, , drop = FALSE]
  ggplot(d, aes(x = gene_order, y = contribution, fill = bar_sign)) +
    geom_col(width = 0.82, color = NA, alpha = 0.96) +
    geom_hline(yintercept = 0, color = "#BFC5D0", linewidth = 0.42) +
    geom_text(
      data = h,
      aes(x = gene_order, y = contribution, label = header_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.25,
      color = INK,
      fontface = "bold",
      lineheight = 0.90
    ) +
    geom_text(
      data = subset(d, show_label),
      aes(
        y = contribution + ifelse(contribution >= 0, ylim * 0.035, -ylim * 0.035),
        label = gene
      ),
      angle = 90,
      size = 1.82,
      fontface = "bold",
      vjust = ifelse(subset(d, show_label)$contribution >= 0, 0, 1),
      color = ifelse(subset(d, show_label)$contribution >= 0, RED, BLUE_NEG),
      show.legend = FALSE
    ) +
    geom_text(
      data = transform(unique(d[, c("display_label", "model"), drop = FALSE]),
                       gene_order = ifelse(model_name == "marine", 70.3, 100.3),
                       contribution = -ylim * 0.80,
                       n_label = paste0("n=", n_pred)),
      aes(x = gene_order, y = contribution, label = n_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0.5,
      size = 2.1,
      color = MUTED
    ) +
    facet_grid(display_label ~ ., switch = "y") +
    scale_fill_manual(
      values = c("positive contribution" = RED, "negative contribution" = BLUE_NEG),
      breaks = c("positive contribution", "negative contribution")
    ) +
    coord_cartesian(ylim = c(-ylim, ylim), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Genes in fixed full-data coefficient order",
      y = "Contribution to fitted logit",
      fill = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", color = INK, size = 13),
      plot.subtitle = element_text(color = MUTED, size = 9),
      axis.title.x = element_text(color = INK, size = 10),
      axis.title.y = element_text(color = INK, size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = GRID, linewidth = 0.25),
      strip.text.y.left = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.52, "lines"),
      legend.position = "bottom",
      legend.justification = "left",
      legend.text = element_text(size = 8.5, color = INK),
      legend.key.size = unit(0.22, "in"),
      plot.margin = margin(10, 12, 10, 13)
    )
}

marine_plot <- plot_one_axis(
  "marine",
  "B. Marine-specialization genomic fingerprints",
  "all 71 predictors; labels mark strongest/story terms",
  71
)
aquatic_plot <- plot_one_axis(
  "binary_aquatic_dependence",
  "C. Aquatic-dependence genomic fingerprints",
  "all 101 predictors; contribution sign is model-score direction",
  101
)

combined <- arrangeGrob(
  marine_plot,
  aquatic_plot,
  ncol = 2,
  top = textGrob(
    "Supplementary Fig. S3 | Internal-branch genomic fingerprints",
    gp = gpar(fontsize = 15, fontface = "bold", col = INK),
    x = unit(0.01, "npc"),
    just = "left"
  ),
  bottom = textGrob(
    "Internal branches are descriptive projected genomic profiles; bars omit the intercept, while profile p includes the intercept.",
    gp = gpar(fontsize = 8.5, col = MUTED)
  )
)

plot_one_axis_subset <- function(model_name, title, subtitle, n_pred, labels_keep, show_legend = FALSE) {
  d <- fp[fp$model == model_name & as.character(fp$display_label) %in% labels_keep, , drop = FALSE]
  d$display_label <- factor(as.character(d$display_label), levels = labels_keep)
  h <- header_df[header_df$model == model_name & as.character(header_df$display_label) %in% labels_keep, , drop = FALSE]
  h$display_label <- factor(as.character(h$display_label), levels = labels_keep)
  ggplot(d, aes(x = gene_order, y = contribution, fill = bar_sign)) +
    geom_col(width = 0.82, color = NA, alpha = 0.96) +
    geom_hline(yintercept = 0, color = "#BFC5D0", linewidth = 0.42) +
    geom_text(
      data = h,
      aes(x = gene_order, y = ylim * 0.80, label = header_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      color = INK,
      fontface = "bold",
      lineheight = 0.90
    ) +
    geom_text(
      data = subset(d, show_label),
      aes(
        y = contribution + ifelse(contribution >= 0, ylim * 0.045, -ylim * 0.045),
        label = gene
      ),
      angle = 90,
      size = 2.45,
      fontface = "bold",
      vjust = ifelse(subset(d, show_label)$contribution >= 0, 0, 1),
      color = ifelse(subset(d, show_label)$contribution >= 0, RED, BLUE_NEG),
      show.legend = FALSE
    ) +
    geom_text(
      data = transform(unique(d[, c("display_label", "model"), drop = FALSE]),
                       gene_order = ifelse(model_name == "marine", 70.3, 100.3),
                       contribution = -ylim * 0.80,
                       n_label = paste0("n=", n_pred)),
      aes(x = gene_order, y = contribution, label = n_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0.5,
      size = 2.6,
      color = MUTED
    ) +
    facet_grid(display_label ~ ., switch = "y") +
    scale_fill_manual(
      values = c("positive contribution" = RED, "negative contribution" = BLUE_NEG),
      breaks = c("positive contribution", "negative contribution")
    ) +
    coord_cartesian(ylim = c(-ylim, ylim), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Genes in fixed full-data coefficient order",
      y = "Contribution to fitted logit",
      fill = NULL
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", color = INK, size = 13),
      plot.subtitle = element_text(color = MUTED, size = 9.2),
      axis.title.x = element_text(color = INK, size = 10),
      axis.title.y = element_text(color = INK, size = 9),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = GRID, linewidth = 0.25),
      strip.text.y.left = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.82, "lines"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.justification = "left",
      legend.text = element_text(size = 8.5, color = INK),
      legend.key.size = unit(0.22, "in"),
      plot.margin = margin(8, 13, 8, 13)
    )
}

row_groups <- list(
  "Cetacean backbone" = c("Cetacea + hippo ancestor", "Crown Cetacea", "Mysticeti", "Odontoceti"),
  "Pinniped split I" = c("Crown Pinnipedia", "Phocidae", "Otarioidea"),
  "Pinniped split II" = c("Otariidae", "Otariinae"),
  "Sirenian and edge context" = c("Sirenia", "Dugong + Hydrodamalis", "Brown bear + polar bear", "Lutrinae / Enhydra-related")
)

group_grobs <- lapply(seq_along(row_groups), function(i) {
  group_name <- names(row_groups)[i]
  labels_keep <- row_groups[[i]]
  left <- plot_one_axis_subset(
    "marine",
    if (i == 1) "B. Marine-specialization genomic fingerprints" else "",
    if (i == 1) "all 71 predictors; labels mark strongest/story terms" else "",
    71,
    labels_keep,
    show_legend = i == length(row_groups)
  )
  right <- plot_one_axis_subset(
    "binary_aquatic_dependence",
    if (i == 1) "C. Aquatic-dependence genomic fingerprints" else "",
    if (i == 1) "all 101 predictors; contribution sign is model-score direction" else "",
    101,
    labels_keep,
    show_legend = i == length(row_groups)
  )
  arrangeGrob(
    textGrob(group_name, x = unit(0.01, "npc"), just = "left", gp = gpar(fontsize = 11, fontface = "bold", col = MUTED)),
    arrangeGrob(left, right, ncol = 2),
    ncol = 1,
    heights = c(0.06, 0.94)
  )
})

four_row <- arrangeGrob(
  grobs = group_grobs,
  ncol = 1,
  top = textGrob(
    "Supplementary Fig. S3 | Internal-branch genomic fingerprints",
    gp = gpar(fontsize = 15, fontface = "bold", col = INK),
    x = unit(0.01, "npc"),
    just = "left"
  ),
  bottom = textGrob(
    "Internal branches are descriptive projected genomic profiles; bars omit the intercept, while profile p includes the intercept.",
    gp = gpar(fontsize = 8.5, col = MUTED)
  )
)

pdf_out <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_simple_two_column.pdf")
png_out <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_simple_two_column.png")
svg_out <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_simple_two_column.svg")

pdf(pdf_out, width = 15.6, height = 22.0, onefile = TRUE)
grid.newpage()
grid.draw(combined)
dev.off()

png(png_out, width = 15.6, height = 22.0, units = "in", res = 300)
grid.newpage()
grid.draw(combined)
dev.off()

svglite(svg_out, width = 15.6, height = 22.0)
grid.newpage()
grid.draw(combined)
dev.off()

label_out <- fp[fp$show_label, c(
  "display_label", "branch_id", "model", "gene", "gene_order", "module",
  "coefficient", "scaled_GBI", "contribution", "abs_contribution", "label_rank"
), drop = FALSE]
write.csv(label_out, file.path(out_dir, "SourceData_ancestor_simple_two_column_displayed_gene_labels.csv"), row.names = FALSE, quote = TRUE, na = "")

note <- c(
  "# Simple Two-Column Ancestor Fingerprint Draft",
  "",
  "This simplified draft removes the profile-map and module-heatmap panels and keeps only the two Figure 6B/C-style fingerprint columns.",
  "",
  "- Left column: marine-specialization axis, all 71 corrected full-data predictors.",
  "- Right column: binary aquatic-dependence axis, all 101 corrected full-data predictors.",
  "- Rows: selected internal branches from the ancestor target map.",
  "- Gene labels: top seven absolute contributors per internal branch and model.",
  "- Contributions: scaled GBI x corrected full-data coefficient; the intercept is omitted from bars but included in the displayed profile probability.",
  "",
  "Evidence boundary: internal branches are descriptive projected genomic profiles, not direct ancestral habitat assignments."
)
writeLines(note, file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_simple_two_column_note.md"))

pdf_four <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_row.pdf")
png_four <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_row.png")
svg_four <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_row.svg")

pdf(pdf_four, width = 16.4, height = 20.0, onefile = TRUE)
grid.newpage()
grid.draw(four_row)
dev.off()

png(png_four, width = 16.4, height = 20.0, units = "in", res = 300)
grid.newpage()
grid.draw(four_row)
dev.off()

svglite(svg_four, width = 16.4, height = 20.0)
grid.newpage()
grid.draw(four_row)
dev.off()

cat("Wrote:\n", pdf_out, "\n", png_out, "\n", svg_out, "\n", sep = "")
cat("Wrote four-row version:\n", pdf_four, "\n", png_four, "\n", svg_four, "\n", sep = "")

plot_card_axis <- function(model_name, labels_keep, title, subtitle, n_pred, show_y_title = FALSE, show_legend = FALSE) {
  d <- fp[fp$model == model_name & as.character(fp$display_label) %in% labels_keep, , drop = FALSE]
  d$display_label <- factor(as.character(d$display_label), levels = labels_keep)
  h <- header_df[header_df$model == model_name & as.character(header_df$display_label) %in% labels_keep, , drop = FALSE]
  h$display_label <- factor(as.character(h$display_label), levels = labels_keep)
  ggplot(d, aes(x = gene_order, y = contribution, fill = bar_sign)) +
    geom_col(width = 0.82, color = NA, alpha = 0.96) +
    geom_hline(yintercept = 0, color = "#BFC5D0", linewidth = 0.38) +
    geom_text(
      data = h,
      aes(x = gene_order, y = ylim * 0.82, label = header_label),
      inherit.aes = FALSE,
      hjust = 0,
      vjust = 1,
      size = 2.45,
      color = INK,
      fontface = "bold",
      lineheight = 0.88
    ) +
    geom_text(
      data = subset(d, show_label),
      aes(
        y = contribution + ifelse(contribution >= 0, ylim * 0.045, -ylim * 0.045),
        label = gene
      ),
      angle = 90,
      size = 2.0,
      fontface = "bold",
      vjust = ifelse(subset(d, show_label)$contribution >= 0, 0, 1),
      color = ifelse(subset(d, show_label)$contribution >= 0, RED, BLUE_NEG),
      show.legend = FALSE
    ) +
    geom_text(
      data = transform(unique(d[, c("display_label", "model"), drop = FALSE]),
                       gene_order = ifelse(model_name == "marine", 70.3, 100.3),
                       contribution = -ylim * 0.82,
                       n_label = paste0("n=", n_pred)),
      aes(x = gene_order, y = contribution, label = n_label),
      inherit.aes = FALSE,
      hjust = 1,
      vjust = 0.5,
      size = 2.15,
      color = MUTED
    ) +
    facet_grid(display_label ~ ., switch = "y") +
    scale_fill_manual(
      values = c("positive contribution" = RED, "negative contribution" = BLUE_NEG),
      breaks = c("positive contribution", "negative contribution")
    ) +
    coord_cartesian(ylim = c(-ylim, ylim), clip = "off") +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Genes in fixed full-data coefficient order",
      y = if (show_y_title) "Contribution to fitted logit" else NULL,
      fill = NULL
    ) +
    theme_minimal(base_size = 8.7) +
    theme(
      plot.title = element_text(face = "bold", color = INK, size = 10.5),
      plot.subtitle = element_text(color = MUTED, size = 7.2),
      axis.title.x = element_text(color = INK, size = 7.6),
      axis.title.y = element_text(color = INK, size = 8.0),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = GRID, linewidth = 0.23),
      strip.text.y.left = element_blank(),
      strip.background = element_blank(),
      panel.spacing.y = unit(0.44, "lines"),
      legend.position = if (show_legend) "bottom" else "none",
      legend.justification = "left",
      legend.text = element_text(size = 6.5, color = INK),
      legend.key.size = unit(0.16, "in"),
      plot.margin = margin(7, 8, 7, 8)
    )
}

make_card <- function(group_name, labels_keep, show_y_title = FALSE) {
  marine_card <- plot_card_axis(
    "marine",
    labels_keep,
    "B. Marine-specialization genomic fingerprints",
    "all 71 predictors; labels mark strongest/story terms",
    71,
    show_y_title = show_y_title,
    show_legend = TRUE
  )
  aquatic_card <- plot_card_axis(
    "binary_aquatic_dependence",
    labels_keep,
    "C. Aquatic-dependence genomic fingerprints",
    "all 101 predictors; contribution sign is model-score direction",
    101,
    show_y_title = show_y_title,
    show_legend = TRUE
  )
  arrangeGrob(
    textGrob(group_name, x = unit(0.01, "npc"), just = "left", gp = gpar(fontsize = 9, fontface = "bold", col = MUTED)),
    arrangeGrob(marine_card, aquatic_card, ncol = 2),
    ncol = 1,
    heights = c(0.035, 0.965)
  )
}

four_card_grobs <- lapply(seq_along(row_groups), function(i) {
  make_card(names(row_groups)[i], row_groups[[i]], show_y_title = i == 1)
})

four_cards <- arrangeGrob(
  grobs = four_card_grobs,
  ncol = 4,
  top = textGrob(
    "Supplementary Fig. S3 | Internal-branch genomic fingerprints",
    gp = gpar(fontsize = 15, fontface = "bold", col = INK),
    x = unit(0.01, "npc"),
    just = "left"
  ),
  bottom = textGrob(
    "Internal branches are descriptive projected genomic profiles; bars omit the intercept, while profile p includes the intercept.",
    gp = gpar(fontsize = 8.2, col = MUTED)
  )
)

pdf_cards <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_cards.pdf")
png_cards <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_cards.png")
svg_cards <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_cards.svg")

pdf(pdf_cards, width = 22.0, height = 10.0, onefile = TRUE)
grid.newpage()
grid.draw(four_cards)
dev.off()

png(png_cards, width = 22.0, height = 10.0, units = "in", res = 300)
grid.newpage()
grid.draw(four_cards)
dev.off()

svglite(svg_cards, width = 22.0, height = 10.0)
grid.newpage()
grid.draw(four_cards)
dev.off()

cat("Wrote four-card version:\n", pdf_cards, "\n", png_cards, "\n", svg_cards, "\n", sep = "")

column_groups <- list(
  "Cetacean + pinniped backbone" = c(
    "Cetacea + hippo ancestor", "Crown Cetacea", "Mysticeti", "Odontoceti",
    "Crown Pinnipedia", "Phocidae", "Otarioidea"
  ),
  "Pinniped + edge context" = c(
    "Otariidae", "Otariinae",
    "Sirenia", "Dugong + Hydrodamalis", "Brown bear + polar bear", "Lutrinae / Enhydra-related"
  )
)

column_grobs <- list(
  plot_card_axis(
    "marine",
    column_groups[[1]],
    "B. Marine-specialization genomic fingerprints",
    "all 71 predictors; labels mark strongest/story terms",
    71,
    show_y_title = TRUE,
    show_legend = TRUE
  ),
  plot_card_axis(
    "binary_aquatic_dependence",
    column_groups[[1]],
    "C. Aquatic-dependence genomic fingerprints",
    "all 101 predictors; contribution sign is model-score direction",
    101,
    show_y_title = TRUE,
    show_legend = TRUE
  ),
  plot_card_axis(
    "marine",
    column_groups[[2]],
    "B. Marine-specialization genomic fingerprints",
    "all 71 predictors; labels mark strongest/story terms",
    71,
    show_y_title = FALSE,
    show_legend = TRUE
  ),
  plot_card_axis(
    "binary_aquatic_dependence",
    column_groups[[2]],
    "C. Aquatic-dependence genomic fingerprints",
    "all 101 predictors; contribution sign is model-score direction",
    101,
    show_y_title = FALSE,
    show_legend = TRUE
  )
)

four_columns <- arrangeGrob(
  grobs = column_grobs,
  ncol = 4,
  top = textGrob(
    "Supplementary Fig. S3 | Internal-branch genomic fingerprints",
    gp = gpar(fontsize = 15, fontface = "bold", col = INK),
    x = unit(0.01, "npc"),
    just = "left"
  ),
  bottom = textGrob(
    "Internal branches are descriptive projected genomic profiles; bars omit the intercept, while profile p includes the intercept.",
    gp = gpar(fontsize = 8.2, col = MUTED)
  )
)

pdf_cols <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns.pdf")
png_cols <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns.png")
svg_cols <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns.svg")

pdf(pdf_cols, width = 18.0, height = 13.0, onefile = TRUE)
grid.newpage()
grid.draw(four_columns)
dev.off()

png(png_cols, width = 18.0, height = 13.0, units = "in", res = 300)
grid.newpage()
grid.draw(four_columns)
dev.off()

svglite(svg_cols, width = 18.0, height = 13.0)
grid.newpage()
grid.draw(four_columns)
dev.off()

cat("Wrote four-column version:\n", pdf_cols, "\n", png_cols, "\n", svg_cols, "\n", sep = "")

pdf_cols_tall <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns_tall.pdf")
png_cols_tall <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns_tall.png")
svg_cols_tall <- file.path(out_dir, "SupplementaryFig_ancestor_fingerprints_four_columns_tall.svg")

pdf(pdf_cols_tall, width = 16.0, height = 15.0, onefile = TRUE)
grid.newpage()
grid.draw(four_columns)
dev.off()

png(png_cols_tall, width = 16.0, height = 15.0, units = "in", res = 300)
grid.newpage()
grid.draw(four_columns)
dev.off()

svglite(svg_cols_tall, width = 16.0, height = 15.0)
grid.newpage()
grid.draw(four_columns)
dev.off()

cat("Wrote tall four-column version:\n", pdf_cols_tall, "\n", png_cols_tall, "\n", svg_cols_tall, "\n", sep = "")
