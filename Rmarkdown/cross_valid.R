#' A more-or-less general cross-validation framework
#' 
#' @param est_func rslt <- est_func(lambda, locs, X, y, NNarray) returns a result list that can be used by pred_func
#' @param pred_func yhat <- pred_func(rslt, locs_pred, X_pred, locs_obs, X_obs, y_obs) returns a vector of predicted values at locs
#' @param crit_func crit_func(y, yhat) returns the loss of the predicted vector
#' @param locs the n-by-d location matrix 
#' @param NNarray the nearest neighborhood array, #col should be at least 1 + m + max(length(testSet))
#' @param X the n-by-p design matrix 
#' @param y the length-n response vector
#' @param m number of neighbors 
#' @param lambda the penalty parameter(s)
#' @param nfold implement the nfold cross validation
#' @param idxRnd a random vector of length n used for selecting the testing dataset
cross_valid <- function(est_func, pred_func, crit_func, locs, NNarray, X, y, m,
                        lambda, nfold, idxRnd = NULL)
{
  n <- length(y)
  if(is.null(idxRnd))
    idxRnd <- sample(c(1 : n), n, F)
  segLen <- floor(n / nfold)
  testIdx <- list()
  for(i in (1 : (nfold - 1)))
    testIdx[[i]] <- idxRnd[((i - 1) * segLen + 1) : (i * segLen)]
  testIdx[[nfold]] <- idxRnd[(i * segLen + 1) : n]
  critSum <- 0
  for(i in 1 : nfold)
  {
    testSet <- testIdx[[i]]
    trainSet <- setdiff(c(1 : n), testSet)
    NNarrayTrain <- NNarray[trainSet, , drop = F]
    NNarrayTrain <- apply(NNarrayTrain, 1, function(x){x[x %in% trainSet][1 : (m + 1)]})
    stopifnot(is.array(NNarrayTrain))
    NNarrayTrain <- t(NNarrayTrain)
    NNarrayTrain <- matrix(sapply(NNarrayTrain, function(x){if(is.na(x)) NA else which(x == trainSet)}),
                           nrow(NNarrayTrain), ncol(NNarrayTrain))
    rslt <- est_func(lambda, locs[trainSet, , drop = F], X[trainSet, , drop = F], 
                     y[trainSet], NNarrayTrain)
    yhat <- pred_func(rslt, locs[testSet, , drop = F], X[testSet, , drop = F], locs[trainSet, , drop = F], 
                      X[trainSet, , drop = F], y[trainSet])
    critSum <- critSum + crit_func(y[testSet], yhat)
  }
  return(critSum / nfold)
}

crit_MSE <- function(y, yhat)
{
  return(sum((y - yhat)^2) / length(y))
}






