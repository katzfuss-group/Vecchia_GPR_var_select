#' A more-or-less general cross-validation framework
#' 
#' @param est_func rslt <- est_func(lambda, locs, X, y, NNarray) returns a result list that can be used by pred_func
#' @param pred_func yhat <- pred_func(rslt, locs_pred, X_pred, locs_obs, X_obs, y_obs) returns a vector of predicted values at locs
#' @param crit_func crit_func(y, yhat) returns the loss of the predicted vector
#' @param locs the n-by-d location matrix 
#' @param X the n-by-p design matrix 
#' @param y the length-n response vector
#' @param NNarray A matrix of indices, usually the output from find_ordered_nn. Row i contains the indices of the observations that observation i conditions on. By convention, the first element of row i is i
#' @param lambda the penalty parameter(s)
#' @param nfold implement the nfold cross validation
#' @param idxRnd a random vector of length n used for selecting the testing dataset
cross_valid <- function(est_func, pred_func, crit_func, locs, X, y, NNarray,
                        lambda, nfold, idxRnd = NULL)
{
  n <- length(y)
  if(is.null(idxRnd))
    idxRnd <- sample(c(1 : n), n, F)
  segLen <- floor(n / nfold)
  testIdx <- list()
  for(i in (1 : (nfold - 1)))
    testIdx[[i]] <- idxRnd[((i - 1) * segLen + 1) : (i * segLen)]
  testIdx[[nfold]] <- idxRnd[((i - 1) * segLen + 1) : n]
  critSum <- 0
  for(i in 1 : nfold)
  {
    testSet <- testIdx[[i]]
    trainSet <- setdiff(c(1 : n), testSet)
    rslt <- est_func(lambda, locs[trainSet, ], X[trainSet, ], 
                     y[trainSet], NNarray[trainSet, ])
    yhat <- pred_func(rslt, locs[testSet, ], X[testSet, ], locs[trainSet, ], X[trainSet, ], 
                      y[trainSet])
    critSum <- critSum + crit_func(y, yhat)
  }
  return(critSum / nfold)
}

crit_MSE <- function(y, yhat)
{
  return(sum((y - yhat)^2) / length(y))
}






