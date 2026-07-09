#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(pROC)
})

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Unexpected argument: ", key)
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    out[[sub("^--", "", key)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

default_path <- function(...) {
  file.path(Sys.getenv("MARINE_MAMMAL_FIGURE2_ROOT", unset = "."), ...)
}

cli_args <- commandArgs(trailingOnly = TRUE)
args <- if (length(cli_args) == 0) list() else parse_args(cli_args)

# Defaults let the script run directly in RStudio with Source.
roc_file <- args$roc %||% default_path("inputs", "Figure2B_gLOOCV_ROC_input.tsv")
pred_file <- args$pred %||% default_path("inputs", "Figure2C_species_prediction_input.tsv")
outdir <- args$outdir %||% default_path("Figure2BC_output")

# Edit these values interactively in RStudio when tuning figure dimensions.
fig2b_width <- 5.8
fig2b_height <- 4.2
fig2c_width <- 7.2
fig2c_height <- 5.8
output_dpi <- 300

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

check_required <- function(data, required, label) {
  missing <- setdiff(required, colnames(data))
  if (length(missing) > 0) {
    stop(label, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

model_colors <- c(
  "marine" = "#2C7FB8",
  "aquatic" = "#D95F0E"
)

model_labels <- c(
  "marine" = "Marine model",
  "aquatic" = "Aquatic-dependence model"
)

theme_nature <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      legend.title = element_blank(),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(colour = "grey30", hjust = 0),
      plot.margin = margin(8, 10, 8, 8)
    )
}

make_roc_data <- function(x) {
  split_x <- split(x, x$model)
  roc_rows <- list()
  auc_rows <- list()
  for (model in names(split_x)) {
    dat <- split_x[[model]]
    dat <- dat[!is.na(dat$true_label) & !is.na(dat$predicted_prob), , drop = FALSE]
    if (length(unique(dat$true_label)) != 2) {
      stop("ROC model '", model, "' does not contain both 0 and 1 labels.")
    }
    roc_obj <- pROC::roc(
      response = dat$true_label,
      predictor = dat$predicted_prob,
      levels = c(0, 1),
      direction = "<",
      quiet = TRUE
    )
    roc_rows[[model]] <- data.frame(
      model = model,
      fpr = 1 - roc_obj$specificities,
      tpr = roc_obj$sensitivities,
      stringsAsFactors = FALSE
    )
    auc_rows[[model]] <- data.frame(
      model = model,
      auc = as.numeric(pROC::auc(roc_obj)),
      stringsAsFactors = FALSE
    )
  }
  list(
    roc = do.call(rbind, roc_rows),
    auc = do.call(rbind, auc_rows)
  )
}

roc_input <- read.delim(roc_file, stringsAsFactors = FALSE, check.names = FALSE)
check_required(
  roc_input,
  c("model", "true_label", "predicted_prob"),
  "ROC input table"
)
roc_input$model <- factor(roc_input$model, levels = c("marine", "aquatic"))

roc_data <- make_roc_data(roc_input)
auc_labels <- setNames(
  sprintf("%s (AUC = %.3f)", model_labels[roc_data$auc$model], roc_data$auc$auc),
  roc_data$auc$model
)

p_roc <- ggplot(roc_data$roc, aes(x = fpr, y = tpr, colour = model)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey72", linewidth = 0.6) +
  geom_path(linewidth = 1.25, lineend = "round") +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_colour_manual(values = model_colors, labels = auc_labels, drop = FALSE) +
  labs(
    title = "Genus-level leave-one-out performance",
    x = "False positive rate",
    y = "True positive rate"
  ) +
  theme_nature(base_size = 11) +
  theme(
    legend.position = c(0.66, 0.22),
    legend.justification = c(0, 0),
    legend.text = element_text(size = 10)
  )

species_input <- read.delim(pred_file, stringsAsFactors = FALSE, check.names = FALSE)
check_required(
  species_input,
  c("species", "common_name", "group", "marine_prob", "aquatic_prob"),
  "Species prediction table"
)

expected_species <- c(
  "Orcinus_orca",
  "Leptonychotes_weddellii",
  "Dugong_dugon",
  "Ursus_maritimus",
  "Enhydra_lutris_kenyoni",
  "Platanista_minor",
  "Hippopotamus_amphibius",
  "Neovison_vison",
  "Ondatra_zibethicus",
  "Lontra_canadensis",
  "Hydrochoerus_hydrochaeris"
)

missing_species <- setdiff(expected_species, species_input$species)
if (length(missing_species) > 0) {
  warning("Selected species missing from prediction table: ", paste(missing_species, collapse = ", "))
}

species_input <- species_input[species_input$species %in% expected_species, , drop = FALSE]
if ("display_order" %in% colnames(species_input)) {
  species_input <- species_input[order(species_input$display_order), , drop = FALSE]
} else {
  species_input$display_order <- match(species_input$species, expected_species)
  species_input <- species_input[order(species_input$display_order), , drop = FALSE]
}

group_levels <- c("Marine category", "Ecologically complex marine case", "Non-marine aquatic")
species_input$group <- factor(species_input$group, levels = group_levels)
species_input$common_name <- factor(
  species_input$common_name,
  levels = rev(species_input$common_name)
)
species_input$y <- as.numeric(species_input$common_name)

species_long <- rbind(
  data.frame(
    species = species_input$species,
    common_name = species_input$common_name,
    group = species_input$group,
    model = "marine",
    predicted_prob = species_input$marine_prob,
    y = species_input$y,
    stringsAsFactors = FALSE
  ),
  data.frame(
    species = species_input$species,
    common_name = species_input$common_name,
    group = species_input$group,
    model = "aquatic",
    predicted_prob = species_input$aquatic_prob,
    y = species_input$y,
    stringsAsFactors = FALSE
  )
)
species_long$model <- factor(species_long$model, levels = c("marine", "aquatic"))

group_ranges <- aggregate(
  y ~ group,
  data = species_input,
  FUN = function(z) c(ymin = min(z) - 0.5, ymax = max(z) + 0.5)
)
group_ranges <- data.frame(
  group = group_ranges$group,
  ymin = group_ranges$y[, "ymin"],
  ymax = group_ranges$y[, "ymax"],
  stringsAsFactors = FALSE
)
group_ranges$fill <- c("#EAF3FB", "#F5F5F5", "#FFF3E6")[match(group_ranges$group, group_levels)]

p_species <- ggplot() +
  geom_rect(
    data = group_ranges,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax, fill = group),
    alpha = 0.55,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = species_input,
    aes(x = marine_prob, xend = aquatic_prob, y = y, yend = y),
    colour = "grey60",
    linewidth = 0.45
  ) +
  geom_point(
    data = species_long,
    aes(x = predicted_prob, y = y, colour = model),
    size = 3.0
  ) +
  scale_y_continuous(
    breaks = species_input$y,
    labels = as.character(species_input$common_name),
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    labels = c("0", "0.25", "0.50", "0.75", "1.00"),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_colour_manual(
    values = model_colors,
    labels = c("Marine model probability", "Aquatic-dependence model probability"),
    drop = FALSE
  ) +
  scale_fill_manual(values = setNames(group_ranges$fill, group_ranges$group), guide = "none") +
  labs(
    title = "Representative species prediction profiles",
    subtitle = "Each row shows predicted probability under the marine and aquatic-dependence models.",
    x = "Predicted probability",
    y = NULL
  ) +
  theme_nature(base_size = 11) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.x = element_line(colour = "grey88", linewidth = 0.35),
    legend.position = "top",
    legend.justification = "right",
    plot.title = element_text(size = 13, face = "bold")
  )

ggsave(file.path(outdir, "Fig2B_ROC.pdf"), p_roc, width = fig2b_width, height = fig2b_height, units = "in", useDingbats = FALSE)
ggsave(file.path(outdir, "Fig2B_ROC.png"), p_roc, width = fig2b_width, height = fig2b_height, units = "in", dpi = output_dpi)
ggsave(file.path(outdir, "Fig2C_species_prediction.pdf"), p_species, width = fig2c_width, height = fig2c_height, units = "in", useDingbats = FALSE)
ggsave(file.path(outdir, "Fig2C_species_prediction.png"), p_species, width = fig2c_width, height = fig2c_height, units = "in", dpi = output_dpi)

message("Figure 2B/2C outputs written to: ", normalizePath(outdir, mustWork = FALSE))
