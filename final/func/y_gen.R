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