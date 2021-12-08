library(scoringRules)
library(FNN)

#' Make independent predictions at locsOut
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @param locsIn in-sample locs
#' @param locsOut out-sample locs
#' @param yIn
#' @param m conditioning set size in scaled Vecchia prediction
#' @param covFn cov kernel name
pred_idp <- function(theta, locsIn, locsOut, yIn, m, covFn)
{
  nIn <- nrow(locsIn)
  d <- ncol(locsIn)
  nOut <- nrow(locsOut)
  idxLocsRel <- which(theta[2 : (d + 1)] > 0)
  dRel <- length(idxLocsRel)
  locsInRel <- locsIn[, idxLocsRel, drop = F]
  locsOutRel <- locsOut[, idxLocsRel, drop = F]
  thetaRel <- theta[c(1, idxLocsRel + 1, d + 2)]
  locsInRelScal <- locsInRel %*% 
    diag(sqrt(thetaRel[2 : (dRel + 1)]), dRel, dRel)
  locsOutRelScal <- locsOutRel %*% 
    diag(sqrt(thetaRel[2 : (dRel + 1)]), dRel, dRel)
  NNarray <- get.knnx(locsInRelScal, locsOutRelScal, m)$nn.index
  NNarray <- cbind(NNarray, (nIn + 1) : (nIn + nOut))
  locsRel <- rbind(locsInRel, locsOutRel)
  mus <- rep(NA, nOut)
  sds <- rep(NA, nOut)
  for(i in 1 : nOut){
    NN <- NNarray[i, ]
    K <- get(covFn)(thetaRel, locsRel[NN, , drop = F])
    L <- t(chol(K))
    mus[i] <- L[m + 1, 1 : m] %*% 
      forwardsolve(L[1 : m, 1 : m], yIn[NN[1 : m]])
    sds[i] <- L[m + 1, m + 1]
  }
  list(mean = mus, sd = sds)
}

#' CRPS out-of-sample score function
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @param locsIn in-sample locs
#' @param locsOut out-sample locs
#' @param yIn
#' @param yOut
#' @param m conditioning set size in scaled Vecchia prediction
#' @param covFn cov kernel name
OOS_crps <- function(theta, locsIn, locsOut, yIn, yOut, m, 
                     covFn = "matern25_scaledim_sqrelevance")
{
  predObj = pred_idp(theta, locsIn, locsOut, yIn, m, covFn)
  mean(crps_norm(y = yOut, mean = predObj$mean, sd = predObj$sd))
}

#' RMSE out-of-sample score function
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @param locsIn in-sample locs
#' @param locsOut out-sample locs
#' @param yIn
#' @param yOut
#' @param m conditioning set size in scaled Vecchia prediction
#' @param covFn cov kernel name
OOS_rmse <- function(theta, locsIn, locsOut, yIn, yOut, m, 
                     covFn = "matern25_scaledim_sqrelevance")
{
  predObj = pred_idp(theta, locsIn, locsOut, yIn, m, covFn)
  sqrt(mean((yOut - predObj$mean)^2) / var(yOut))
}


#' Split locs and y to p1 : (1 - p1)
#' 
#' @param locs  n X d matrix
#' @param y n vec
#' @param p1 first part proportion
splt_locs_y <- function(locs, y, p1, seed = NULL){
  if(!is.null(seed))
    set.seed(seed)
  n <- length(y)
  idx1 <- sample(1 : n, round(p1 * n))
  idx2 <- setdiff(1 : n, idx1)
  list(locs1 = locs[idx1, , drop = F], y1 = y[idx1], 
       locs2 = locs[idx2, , drop = F], y2 = y[idx2])
}

#' Reconstruct theta from the positive part of theta
#' 
#' @param thetaRel  positive part of theta
#' @param idx positive indices
#' @param d total # covariates
get_theta <- function(thetaRel, idx, d)
{
  theta <- rep(0, d + 2)
  theta[c(1, idx + 1, d + 2)] <- thetaRel
  theta
}

#' Wrapper for vecchia_meanzero_loglik_grad_info
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @param locs  n X d matrix
#' @param y n vec
#' @param NNarray NNarray
#' @param miniGrad mini-batch size for gradient computation
#' @param covFn covariance function name
comp_grad_mean0 <- function(theta, locs, y, NNarray, miniGrad, 
                      covFn = "matern25_scaledim_sqrelevance")
{
  n <- nrow(NNarray)
  batchIdx <- sample(x = 1 : n, size = miniGrad, replace = F)
  vecchia_meanzero_loglik_grad_info(theta, covFn, y, locs, 
                                    NNarray[batchIdx, , drop = F])
}

#' Wrapper for maximin ordering and NNarray computation
#' 
#' @param theta covariance parms (var, R/SR/range(s), tau)
#' @param idx effective covariate indices
#' @param locs  n X d matrix
#' @param y n vec
#' @param lb lower bounds for theta
#' @param m conditioning set size in scaled Vecchia prediction
MM_NN <- function(theta, idx, locs, y, lb, m)
{
  n <- nrow(locs)
  d <- ncol(locs)
  dHat <- length(idx)
  locsRel <- locs[, idx, drop = F]
  locsRelScal <- locsRel %*% diag(sqrt(theta[1 + idx]), dHat, dHat)
  odr <- order_maxmin(locsRelScal)
  yOdr <- y[odr]
  locsRelScalOdr <- locsRelScal[odr, , drop = F]
  locsOdr <- locs[odr, , drop = F]
  locsOdrRel <- locsOdr[, idx, drop = F] 
  NNarray <- GpGp::find_ordered_nn(locsRelScalOdr, m = m)
  thetaRel <- theta[c(1, idx + 1, d + 2)]
  lbRel <- lb[c(1, idx + 1, d + 2)]
  list(locsOdr = locsOdr, locsOdrRel = locsOdrRel, yOdr = yOdr, 
       thetaRel = thetaRel, lbRel = lbRel, NNarray = NNarray)
}