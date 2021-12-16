

#' Bridge penalty
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @lambda penalty scalar
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

#' Averaged bridge penalty
#' 
#' @param thetaAvg averaged covariance parms (var, R/SR/range(s), tau)
#' @lambda penalty scalar
#' @iter `#iter to average over
pen_avg_brdg <- function(thetaAvg, lambda, iter){
  lambda * sum(thetaAvg[-c(1, length(thetaAvg))]^(0.25))
}
dpen_avg_brdg <- function(thetaAvg, lambda, iter){
  rAvg <- thetaAvg[-c(1, length(thetaAvg))]
  rpen <- lambda * rAvg^(0.25 - 1) * 0.25 / iter
  rpen[rAvg < 1e-20] <- lambda * (1e-20)^(0.25 - 1) * 0.25 / iter
  c(0, rpen, 0)
}
ddpen_avg_brdg <- function(thetaAvg, lambda, iter){
  diag(rep(0, length(thetaAvg)))
}



