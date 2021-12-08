

#' Bridge penalty
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @lambda 
pen_brdg <- function(theta, lambda){
  lambda * sum(theta[-c(1, length(theta))]^(0.25))
}
dpen_brdg <- function(theta, lambda){
  r <- theta[-c(1, length(theta))]
  rpen <- lambda * r^(0.25 - 1) * 0.25
  rpen[r < 1e-20] <- lambda * (1e-20)^(0.25 - 1) * 0.25
  c(0, rpen, 0)
}
ddpen_brdg <- function(theta, lambda){
  diag(rep(0, length(theta)))
}