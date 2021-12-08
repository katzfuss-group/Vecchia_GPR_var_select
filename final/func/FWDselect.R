library(parallel)
source("fisher_scoring.R")
source("misc.R")

model_fit <- function(locsIdx, m, theta, locsTrain, yTrain, covFn, OOS_score, 
                      link_func, resp_func, dresp_func)
{
  maxIter <- 100
  # MM and NN
  MMNNObj <- MM_NN(theta, locsIdx, locsTrain, yTrain, 
                   rep(-Inf, length(theta)), m)
  # locsFit <- locsTrain[, locsIdx, drop = F]
  # thetaFit <- theta[c(1, locsIdx + 1, d + 2)]
  # locsFitScal <- locsFit %*% 
  #   diag(sqrt(resp_func(thetaFit[2 : (1 + length(locsIdx))])), 
  #        length(locsIdx), length(locsIdx))
  # odr <- 1 : nTrain # order_maxmin(locsFitScal)
  # yOdr <- yTrain[odr]
  # locsFitScalOdr <- locsFitScal[odr, , drop = F]
  # locsFitOdr <- locsFit[odr, , drop = F]
  # NNarray <- GpGp::find_ordered_nn(locsFitScalOdr, m = m)
  
  objfun <- function(thetaTrans){
    likobj <- 
      vecchia_meanzero_loglik_grad_info(resp_func(thetaTrans), 
                                        covFn,
                                        MMNNObj$yOdr, MMNNObj$locsOdrRel, 
                                        MMNNObj$NNarray)
    likobj$loglik <- -likobj$loglik
    likobj$grad <- -c(likobj$grad) * dresp_func(thetaTrans)
    likobj$info <- likobj$info * outer(dresp_func(thetaTrans),
                                       dresp_func(thetaTrans))
    return(likobj)
  }
  
  FisherObj <- fisher_scoring(objfun, MMNNObj$thetaRel, resp_func, T, 1e-4, 
                              maxIter)
  thetaFit <- FisherObj$logparms
  thetaTmp <- rep(0, d + 2)
  thetaTmp[c(1, locsIdx + 1, d + 2)] <- resp_func(thetaFit)
  OOSScore <- OOS_score(theta = thetaTmp, m = m)
  list(theta = thetaFit, score = OOSScore)
}

### TBD
forward_selection <- function(m = 100, outputFn){
  idxSet <- list()
  thetaSet <- list()
  scoreSet <- list()
  selectIdx <- c()
  
  theta <- link_func(c(1, rep(1.0, d), 0.0001))
  
  link_func <- log
  resp_func <- exp
  dresp_func <- exp
  
  startTime <- Sys.time()
  for(i in 1 : d)
  {
    cat("\n====================================\n")
    unselectIdx <- setdiff(1 : d, selectIdx)
    idxList <- lapply(unselectIdx, function(x){c(selectIdx, x)})
    models <- mclapply(idxList, model_fit, mc.cores = 8, 
                       m = m, theta = theta)
    # models <- lapply(idxList, model_fit, m = m, theta = theta)
    scores <- unlist(lapply(models, function(x){x$score}))
    idxBest <- which.min(scores)
    cat("i =", i, "idxBest is", idxList[[idxBest]], "best score is", 
        min(scores), "\n")
    # Store results
    idxSet[[i]] <- idxList[[idxBest]] + 1
    thetaSet[[i]] <- rep(0, d + 2)
    thetaSet[[i]][c(1, idxSet[[i]], d + 2)] <- 
      resp_func(models[[idxBest]]$theta)
    scoreSet[[i]] <- models[[idxBest]]$score
    theta[c(1, idxSet[[i]], d + 2)] <- models[[idxBest]]$theta
    selectIdx <- idxList[[idxBest]]
    if(i > 1 && scoreSet[[i]] / scoreSet[[i - 1]] > 0.99 &&
       length(setdiff(idxSet[[i]], idxSet[[i - 1]])) > 0 && 
       sum(thetaSet[[i]][2 : (d + 1)]) > 0.1)
      break
  }
  endTime <- Sys.time()
  timeObj <- endTime - startTime
  timeUsed <- as.numeric(timeObj, units = "secs")
  cat("Forward selection used", timeUsed, "seconds\n")
  save(list = c("idxSet", "thetaSet", "scoreSet", "timeUsed"), 
       file = outputFn)
}