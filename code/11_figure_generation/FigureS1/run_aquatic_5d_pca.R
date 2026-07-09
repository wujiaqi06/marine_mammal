#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: run_aquatic_5d_pca.R pca_input.tsv slim_table.tsv out_dir", call. = FALSE)
}

pca_input <- args[[1]]
slim_table <- args[[2]]
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pca_df <- read.delim(pca_input, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(slim_table, stringsAsFactors = FALSE, check.names = FALSE)

score_cols <- c(
  "score_foraging_medium_0_3",
  "score_locomotion_escape_0_3",
  "score_reproduction_nursery_0_3",
  "score_morphophysiology_0_5",
  "score_time_budget_0_4"
)

scores_matrix <- as.matrix(pca_df[, score_cols])
rownames(scores_matrix) <- pca_df$species_id

pca <- prcomp(scores_matrix, center = TRUE, scale. = TRUE)
variance <- data.frame(
  PC = paste0("PC", seq_along(pca$sdev)),
  eigenvalue = pca$sdev^2,
  variance_explained = (pca$sdev^2) / sum(pca$sdev^2),
  cumulative_variance = cumsum((pca$sdev^2) / sum(pca$sdev^2))
)

loadings <- data.frame(
  trait_dimension = rownames(pca$rotation),
  pca$rotation,
  check.names = FALSE
)

pc_scores <- data.frame(
  species_id = rownames(pca$x),
  pca$x,
  check.names = FALSE
)
pc_scores <- merge(
  pc_scores,
  meta[, c("species_id", "common_name_paper", "aquaticity_score_sum_0_18", "score_status", "score_confidence")],
  by = "species_id",
  all.x = TRUE,
  sort = FALSE
)

write.table(variance, file.path(out_dir, "aquaticity_5d_pca_variance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(loadings, file.path(out_dir, "aquaticity_5d_pca_loadings.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(pc_scores, file.path(out_dir, "aquaticity_5d_pca_scores.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

status_levels <- sort(unique(pc_scores$score_status))
palette <- c(
  clear_terrestrial_zero = "#9E9E9E",
  manual_curated_terrestrial = "#6E6E6E",
  reviewed_r1_initial = "#D55E00",
  marine_obligate_draft = "#0072B2",
  manual_curated_extinct_marine = "#004C7F",
  pinniped_draft = "#009E73",
  marine_semi_aquatic_draft = "#CC79A7"
)
point_cols <- palette[pc_scores$score_status]
point_cols[is.na(point_cols)] <- "#333333"

plot_pca <- function(filename, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    pdf(filename, width = 7.5, height = 6)
  } else {
    png(filename, width = 1800, height = 1450, res = 220)
  }
  par(mar = c(5, 5, 3, 9), xpd = TRUE)
  plot(
    pc_scores$PC1, pc_scores$PC2,
    pch = 21,
    bg = point_cols,
    col = "white",
    cex = 1.15,
    xlab = sprintf("PC1 (%.1f%%)", 100 * variance$variance_explained[1]),
    ylab = sprintf("PC2 (%.1f%%)", 100 * variance$variance_explained[2]),
    main = "Aquaticity 5D PCA"
  )
  abline(h = 0, v = 0, col = "gray85", lty = 2)
  legend(
    "right",
    inset = c(-0.42, 0),
    legend = names(palette)[names(palette) %in% status_levels],
    pt.bg = palette[names(palette) %in% status_levels],
    pch = 21,
    pt.cex = 1.2,
    bty = "n",
    cex = 0.75
  )

  label_species <- c(
    "Hydrodamalis_gigas",
    "Trichechus_manatus_latirostris",
    "Ondatra_zibethicus",
    "Lontra_canadensis",
    "Pteronura_brasiliensis",
    "Bubalus_bubalis",
    "Syncerus_caffer",
    "Alces_alces",
    "Neovison_vison",
    "Ursus_maritimus"
  )
  idx <- which(pc_scores$species_id %in% label_species)
  text(
    pc_scores$PC1[idx],
    pc_scores$PC2[idx],
    labels = pc_scores$common_name_paper[idx],
    pos = 4,
    cex = 0.62,
    offset = 0.25
  )
  dev.off()
}

plot_pca(file.path(out_dir, "aquaticity_5d_pca_PC1_PC2.pdf"), "pdf")
plot_pca(file.path(out_dir, "aquaticity_5d_pca_PC1_PC2.png"), "png")

cat("Wrote PCA outputs to", out_dir, "\n")
