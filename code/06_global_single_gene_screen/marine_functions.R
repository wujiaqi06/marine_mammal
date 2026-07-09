median_no_na <- function(x) {
  x0 <- x[!is.na(x)]
  if (length(x0) == 0) {
    return(NA_real_)
  }
  median(x0)
}

mean_no_na <- function(x) {
  x0 <- x[!is.na(x)]
  if (length(x0) == 0) {
    return(NA_real_)
  }
  mean(x0)
}

trim_upper_975 <- function(x) {
  if (all(is.na(x))) {
    return(x)
  }
  cutoff <- as.numeric(stats::quantile(x, 0.975, na.rm = TRUE, names = FALSE))
  x[x >= cutoff] <- NA
  x
}

apply_margin_preserve_matrix <- function(m, margin, fun) {
  if (!is.matrix(m) && !is.data.frame(m)) {
    stop("Input must be a matrix or data.frame.")
  }
  m <- as.matrix(m)
  if (margin == 1) {
    out <- t(vapply(seq_len(nrow(m)), function(i) fun(m[i, ]), FUN.VALUE = numeric(ncol(m))))
    rownames(out) <- rownames(m)
    colnames(out) <- colnames(m)
    return(out)
  }
  if (margin == 2) {
    out <- vapply(seq_len(ncol(m)), function(i) fun(m[, i]), FUN.VALUE = numeric(nrow(m)))
    rownames(out) <- rownames(m)
    colnames(out) <- colnames(m)
    return(out)
  }
  stop("margin must be 1 or 2.")
}

impute_by_column_mean <- function(matrix0) {
  matrix0 <- as.matrix(matrix0)
  out <- matrix0
  for (i in seq_len(ncol(out))) {
    replace_value <- mean(out[, i], na.rm = TRUE)
    out[is.na(out[, i]), i] <- replace_value
  }
  as.data.frame(out)
}

impute_by_column_terminal_mean <- function(matrix0, terminal_mean) {
  matrix0 <- as.matrix(matrix0)
  out <- matrix0
  for (i in seq_len(ncol(out))) {
    gene0 <- colnames(out)[i]
    if (!gene0 %in% names(terminal_mean)) {
      stop(sprintf("Missing terminal mean for gene: %s", gene0))
    }
    out[is.na(out[, i]), i] <- terminal_mean[[gene0]]
  }
  as.data.frame(out)
}

assert_unique_ids <- function(ids, label) {
  dup <- unique(ids[duplicated(ids)])
  if (length(dup) > 0) {
    stop(sprintf("Duplicated %s detected: %s", label, paste(dup, collapse = ", ")))
  }
}
