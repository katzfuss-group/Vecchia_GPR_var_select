library(mvtnorm)
library(lhs)

#' Generate n X d independent locs
#' 
#' @param n # locs
#' @param d # covariates
locs_gen_idp <- function(n, d, nmlz = T){
  locs <- randomLHS(n, d)
  if(nmlz)
    locs * outer(rep(sqrt(n), n), 
                 1 / sqrt(colSums(locs^2)))
  else
    locs
}

#' Generate n X d dependent locs
#' 
#' @param n # locs
#' @param d # covariates
#' @param rho # corr btw covariates
locs_gen_dp <- function(n, d, rho, nmlz = T){
  covM <- matrix(rho, d, d)
  diag(covM) <- 1
  locs <- rmvnorm(n, sigma = covM) 
  if(nmlz)
    locs * outer(rep(sqrt(n), n), 
                 1 / sqrt(colSums(locs^2)))
  else
    locs
}