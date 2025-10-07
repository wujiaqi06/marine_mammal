
median_noNA <- function(x){
  x0 <- x[!is.na(x)]
  median(x0)
}

mean_noNA <- function(x){
  x0 <- x[!is.na(x)]
  mean(x0)
}

outlier0 <- function(x){
  up_line <- quantile(x, 0.975, na.rm = T)
  out <- x[x >= up_line]
  x[x %in% out] <- NA
  x
}

host_apply <- function(m, x, f){
  m_out <- data.frame(matrix(0,nrow = nrow(m), ncol = ncol(m)))
  rownames(m_out) <- rownames(m)
  colnames(m_out) <- colnames(m)
  if (x == 1){
    for (i in 1:nrow(m)){
      f1 <- rownames(m)[i]
      d0 <- m[i,]
      d1 <- f(d0)
      m_out[f1,] <- d1
    }
  }else if (x == 2){
    for (i in 1:ncol(m)){
      f1 <- colnames(m)[i]
      d0 <- m[,i]
      d1 <- f(d0)
      m_out[,f1] <- d1
    }
  }
  m_out
}

imputation_by_column <- function(matrix0){
  matrix_imputate <- c()
  for (i in 1:ncol(matrix0)){
    imputed_mean = replace(matrix0[,i], is.na(matrix0[,i]), mean(matrix0[,i], na.rm = TRUE))
    matrix_imputate <- rbind(matrix_imputate, imputed_mean)
  }
  dim(matrix_imputate);matrix_imputate[c(1:3),c(1:4)]
  rownames(matrix_imputate) <- colnames(matrix0)
  colnames(matrix_imputate) <- rownames(matrix0)
  matrix_imputate <- t(matrix_imputate)
  data.frame(matrix_imputate)
}

imputation_by_column_terminal <- function(matrix0, terminal_mean){
  matrix_imputate <- c()
  for (i in 1:ncol(matrix0)){
    gene0 <- colnames(matrix0)[i]
    imputed_mean = replace(matrix0[,i], is.na(matrix0[,i]), terminal_mean[gene0])
    matrix_imputate <- rbind(matrix_imputate, imputed_mean)
  }
  dim(matrix_imputate);matrix_imputate[c(1:3),c(1:4)]
  rownames(matrix_imputate) <- colnames(matrix0)
  colnames(matrix_imputate) <- rownames(matrix0)
  matrix_imputate <- t(matrix_imputate)
  data.frame(matrix_imputate)
}
