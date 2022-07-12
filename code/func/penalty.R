

#' Bridge penalty
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @lambda penalty scalar
pen_brdg <- function(theta, lambda, gamma = 0.25){
  lambda * sum(theta[-c(1, length(theta))]^(gamma))
}
dpen_brdg <- function(theta, lambda, gamma = 0.25){
  r <- theta[-c(1, length(theta))]
  rpen <- lambda * r^(gamma - 1) * gamma
  rpen[r < 1e-20] <- lambda * (1e-20)^(gamma - 1) * gamma
  c(0, rpen, 0)
}
ddpen_brdg <- function(theta, lambda, gamma = 0.25){
  diag(rep(0, length(theta)))
}

#' Averaged bridge penalty
#' 
#' @param thetaAvg averaged covariance parms (var, R/SR/range(s), tau)
#' @lambda penalty scalar
#' @iter #iter to average over
pen_avg_brdg <- function(thetaAvg, lambda, iter, gamma = 0.25){
  lambda * sum(thetaAvg[-c(1, length(thetaAvg))]^(gamma))
}
dpen_avg_brdg <- function(thetaAvg, lambda, iter, gamma = 0.25){
  rAvg <- thetaAvg[-c(1, length(thetaAvg))]
  rpen <- lambda * rAvg^(gamma - 1) * gamma / iter
  rpen[rAvg < 1e-20] <- lambda * (1e-20)^(gamma - 1) * gamma / iter
  c(0, rpen, 0)
}
ddpen_avg_brdg <- function(thetaAvg, lambda, iter, gamma = 0.25){
  diag(rep(0, length(thetaAvg)))
}



