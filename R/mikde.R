#' KDE estimator of mutual information
#' @export

mikde <- function(x, y, kden= 100){
  xy <- MASS::kde2d(x = x, y = y, n = kden)
  pxy <- xy$z
  pxy <- (pxy - min(pxy)) / (max(pxy) - min(pxy))
  pxy <- pxy / sum(pxy)

  px <- replicate(ncol(pxy), rowSums(pxy))
  py <- t(replicate(nrow(pxy), colSums(pxy)))

  mi <- sum(pxy * log(pxy / (px * py)), na.rm = T)
  return(mi)
}
