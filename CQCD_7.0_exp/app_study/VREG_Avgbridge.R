penfun <- function(thetaAvg, lambda, iter){
  lambda * sum(thetaAvg[-c(1, length(thetaAvg))]^(0.25))
}
dpenfun <- function(thetaAvg, lambda, iter){
  rAvg <- thetaAvg[-c(1, length(thetaAvg))]
  rpen <- lambda * rAvg^(0.25 - 1) * 0.25 / iter
  rpen[rAvg < 1e-20] <- lambda * (1e-20)^(0.25 - 1) * 0.25 / iter
  c(0, rpen, 0)
}
ddpenfun <- function(thetaAvg, lambda, iter){
  diag(rep(0, length(thetaAvg)))
}

forward_backward <- function(lambda, theta0, idx0, gradObj0, m, k, lb, 
                             arg_check_SR, batchSzCQCD, batchSzGrad)
{
  theta <- theta0
  idx <- idx0
  gradObj <- gradObj0
  OOS <- Inf
  while(T)
  {
    odrDec <- order(gradObj$grad[2 : (d + 1)], decreasing = T)
    if(length(idx) == 0)
      idxNew <- 1 + odrDec[1 : k]
    else
      idxNew <- c(idx, 1 + setdiff(odrDec, idx - 1)[1 : k])
    idxSel <- setdiff(idxNew, idx)
    cat("Selected var:", idxSel, "\n")
    # opt with SR parameters
    thetaRel <- theta[c(1, idxNew, 2 + d)] 
    thetaRel[(length(idx) + 2) : (length(idx) + k + 1)] <- 0.01
    # maximin order and NNarray
    locsRel <- locsTrain[, idxNew - 1, drop = F]
    locsRelScal <- locsRel %*% 
      diag(sqrt(thetaRel[2 : (1 + length(idxNew))]), 
           length(idxNew), length(idxNew))
    odr <- 1 : nTrain # GpGp::order_maxmin(locsRelScal)
    yOdr <- yTrain[odr]
    locsRelScalOdr <- locsRelScal[odr, , drop = F]
    locsOdr <- locsTrain[odr, , drop = F]
    NNarray <- GpGp::find_ordered_nn(locsRelScalOdr, m = m)
    # extract relevant col in locsOdr and lb
    locsOdrRel <- locsOdr[, idxNew - 1, drop = F] 
    lbRel <- lb[c(1, idxNew, 2 + d)] 
    objfun <- function(theta, batchIdx, thetaAvg, iter){
      likObj <- GpGp::vecchia_meanzero_loglik(theta, 
                                              "matern25_scaledim_sqrelevance",
                                              yOdr, locsOdrRel, 
                                              NNarray[batchIdx, , drop = F])
      likObj$loglik <- likObj$loglik - penfun(thetaAvg, lambda, iter)
      likObj
    }
    objfun_gdfm <- function(theta, batchIdx, thetaAvg, iter){
      likObj <- GpGp::vecchia_meanzero_loglik_grad_info(theta, 
                                              "matern25_scaledim_sqrelevance",
                                              yOdr, locsOdrRel, 
                                              NNarray[batchIdx, , drop = F])
      likObj$loglik <- likObj$loglik - penfun(thetaAvg, lambda, iter)
      likObj$grad <- likObj$grad - dpenfun(thetaAvg, lambda, iter)
      likObj$info <- likObj$info + ddpenfun(thetaAvg, lambda, iter)
      likObj
    }
    timeObj <- system.time(
      optObj <- CQCD_stochastic_taper(objfun, objfun_gdfm, nTrain, batchSzCQCD,
                                      thetaRel, maxIterOut = 500, 
                                      maxIterIn = 40, lb = lbRel,
                                      arg_check = arg_check_SR))
    # Take out zero relevance
    thetaRel <- optObj$covparms
    idxLocZero <- which(thetaRel[2 : (1 + length(idxNew))] == 0)
    if(length(idxLocZero) > 0)
    {
      idxDesel <- idxNew[idxLocZero]
      cat("Predictor", idxDesel, "are zerod out \n")
      idxNew <- idxNew[-idxLocZero]
      thetaRel <- thetaRel[-(idxLocZero + 1)]
    }
    thetaNew <- rep(0, length(theta))
    thetaNew[c(1, idxNew, length(theta))] <- thetaRel
    OOSNew <- OOS_score(thetaNew, m)
    if(OOSNew > OOS * 0.99 || length(idx) > 100)
      break
    idx <- idxNew
    theta <- thetaNew
    OOS <- OOSNew
    # compute grad with thetaNew
    batchIdx <- sample(x = 1 : nTrain, size = batchSzGrad, replace = F)
    gradObj <- 
      GpGp::vecchia_meanzero_loglik_grad_info(theta, 
                                              "matern25_scaledim_sqrelevance",
                                              yOdr, locsOdr, 
                                              NNarray[batchIdx, , drop = F])
  }
  cat("lambda =", lambda, "idx =", idx, "\n\n")
  return(list(lambda = lambda, theta = theta, idx = idx, gradObj = gradObj,
              OOS = OOS))
}

VReg <- function(m = 100, k = 5, outputFn = "VREG.RData", 
                 lambdaUp = 16, lambdaLo = 1/128){
  lambdaVec <- c(rev(2^(seq(from = log2(lambdaLo), 
                            to = log2(max(lambdaLo, lambdaUp)), 
                            by = 1))))
  niter <- length(lambdaVec)
  batchSzCQCD <- min(128, nTrain)
  batchSzGrad <- min(128, nTrain)
  
  theta <- c(1, rep(1e-8, d), 0.0001)
  lb <- c(0.01^2, rep(0, d), 0.01^2) 
  arg_check_SR <- function(x) {sum(sqrt(x[-c(1, length(x))])) > 1e-4}
  idx <- c()
  yOdr <- yTrain
  locsOdr <- locsTrain
  NNarray <- find_ordered_nn(locsTrain[, sample(1 : d, k), drop = F], m = m)
  batchIdx <- sample(x = 1 : nTrain, size = batchSzGrad, replace = F)
  gradObj <- vecchia_meanzero_loglik_grad_info(theta, 
                                               "matern25_scaledim_sqrelevance",
                                               yOdr, locsOdr, 
                                               NNarray[batchIdx, , drop = F])
  lb[which.max(gradObj$grad[2 : (d + 1)]) + 1] <- 2e-4
  
  idxSet <- list()
  thetaSet <- list()
  scoreSet <- list()
  lambdaSet <- list()
  loglikSet <- list()
  startTime <- Sys.time()
  for(i in 1 : niter)
  {
    cat("\n====================================\n")
    lambda <- lambdaVec[i]
    cat("i =", i, "lambda =", lambda, "\n")
    # fit model with penalty
    optObj <- forward_backward(lambda, theta, idx, gradObj, m, k, lb, 
                               arg_check_SR, batchSzCQCD, batchSzGrad)
    theta <- optObj$theta
    idx <- optObj$idx
    gradObj <- optObj$gradObj
    # Store results
    idxSet[[i]] <- optObj$idx
    thetaSet[[i]] <- optObj$theta
    scoreSet[[i]] <- optObj$OOS # OOS_score(optObj$theta, m)
    lambdaSet[[i]] <- lambda
    loglikSet[[i]] <- optObj$gradObj$loglik
    cat("OOS score is", scoreSet[[i]], "\n")
    if(i > 1 && scoreSet[[i]] / scoreSet[[i - 1]] > 0.99 &&
       length(setdiff(idxSet[[i]], idxSet[[i - 1]])) > 0 && 
       sum(thetaSet[[i]][2 : (d + 1)]) > 0.1)
      break
    # Reset the variance and nugget parameter
    theta[1] <- 1
    theta[length(theta)] <- 0.0001
  }
  endTime <- Sys.time()
  timeObj <- endTime - startTime
  timeUsed <- as.numeric(timeObj, units = "secs")
  cat("VREG used", timeUsed, "seconds\n")
  save(list = c("idxSet", "thetaSet", "scoreSet", "lambdaSet", "timeUsed"), 
       file = outputFn)
}

