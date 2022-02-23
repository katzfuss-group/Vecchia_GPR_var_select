library(mvtnorm)
library(GpGp)
library(parallel)
library(mvnfast)
source("piston_func.R")

#' Generate y based on the Piston function
#' 
#' @param locs n X d matrix
Pist_gen <- function(locs){
  n <- nrow(locs)
  d <- ncol(locs)
  if(d < 7){
    stop("Piston function takes 7 covariates at least\n")
  }
  locsTmp <- apply(locs[, 1 : 7], 2, function(x){x - min(x) + 1})
  apply(locsTmp, 1, pistonfun)
}

#' Generate multivariate normal y
#' 
#' @param locs n X d matrix
#' @param parms covariance parameters
#' @param covFn covariance function name from 'GpGp'
MVN_gen <- function(locs, parms, covFn){
  n <- nrow(locs)
  covM <- get(covFn)(parms, locs)
  # as.vector(rmvnorm(n = 1, sigma = covM))
  cat("Generating y\n")
  timeUsed <- system.time(
    result <- as.vector(rmvn(n = 1, mu = rep(0, n), 
                             sigma = covM,
                             ncores = detectCores() - 1)))[[3]]
  cat("Generating y used", timeUsed, "seconds\n")
  result
}
