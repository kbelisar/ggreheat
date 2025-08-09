### FROM PThelper https://github.com/PBCAR/PThelper

#' `upper_tri()`
#'
#' @param x An R object.
#' @noRd

upper_tri <- function(x){
  x[lower.tri(x)] <- NA
  return(x)
}

#' `base_melt()`
#'
#' @param x An R object.
#' @noRd

base_melt <- function(data){

  dat_name <- c(colnames(data))

  dat <- data.frame(expand.grid(lapply(dim(data), seq_len)),value = as.vector(as.matrix(data)))
  dat <- dat[!is.na(dat$value),]
  dat$Var1 <- factor(dat$Var1, levels = c(seq(dat_name)), labels = c(dat_name))
  dat$Var2 <- factor(dat$Var2, levels = c(seq(dat_name)), labels = c(dat_name))

  return(dat)
}
