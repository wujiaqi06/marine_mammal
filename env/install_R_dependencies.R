#!/usr/bin/env Rscript

required <- c(
  "ape", "castor", "glmnet", "ggplot2", "gridExtra", "svglite",
  "pROC", "cowplot", "ggrepel", "digest"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}
cat("R dependencies available:", paste(required, collapse = ", "), "\n")
