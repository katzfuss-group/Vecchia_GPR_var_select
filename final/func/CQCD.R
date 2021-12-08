#' Coordinate descent using the 2nd order approximation
#' 
#' @param likfun likelihood function, returns log-likelihood
#' @param likfun_GDFIM returns log-likelihood, gradient, and FIM
#' @param parms0 starting values of parameters
#' @param convtolOut convergence tolerance on step dot grad
#' @param convtolIn convergence tolerance on the step of one coordinate descent epoch
#' @param cAmij the Armijo condition constant
#' @param maxIterOut maximum number of 2nd order approximations
#' @param maxIterIn maximum number of epochs in coordinate descent
#' @param lb the lower bounds for all parameters
#' @param arg_check the function for checking parms
#' @param silent TRUE/FALSE for suppressing output
CQCD <- function(likfun, likfun_GDFIM, parms0, 
                 convtolOut = 1e-4, convtolIn = 1e-4, 
                 cAmij = 1e-4, maxIterOut = 20, maxIterIn = 40, 
                 lb = rep(0, length(parms0)),
                 arg_check = function(...){T}, silent = FALSE)
{
  if(any(parms0 < lb))
    stop(paste("Initial values for parameters ", which(parms0 < lb), 
               "are smaller than their lower bounds. Stopping...  \n"))
  if(!arg_check(parms0))
    stop(paste("Argument check for parms0 failed. Stopping... \n"))
  
  parms <- parms0
  likobj <- likfun_GDFIM(parms)
  # record first loglik value
  loglikInit <- likobj$loglik
  
  for(i in 1 : maxIterOut)
  {
    # check for NAs and NaNs
    if(anyNA(likobj)){
      warning(paste("NA found in likobj at iteration", i - 1, "early return...\n"))
      return(list(covparms = parms, loglik = likobj$loglik, grad = likobj$grad,
                  info = likobj$info, neval = i))
    }
    
    obj <- - likobj$loglik
    grad <- - likobj$grad 
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("CQCD iter ", i - 1, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("\n")
    }
    # coordinate descent
    b <- grad - as.vector(H %*% parms)
    DSCObj <- dsc_lb(H, b, parms, silent, convtolIn, maxIterIn, lb)
    parmsNew <- bktk_Armijo(parms, obj, grad, DSCObj$parms, cAmij, 
                            function(x){- likfun(x)$loglik}, lb, arg_check)
    # if coord descent meets Armijo condition
    if((abs(sum((parmsNew - parms) * grad)) >= convtolOut))
    {
      if(!silent)
        cat("Lower bounded coordinate descent with FIM is used \n")
    }
    else # Use lower bounded gradient descent
    {
      parmsNew <- gd_Armijo(parms, obj, grad, cAmij, 
                            function(x){- likfun(x)$loglik}, lb, arg_check)
      if(!silent)
        cat("Lower bounded gradient descent is used \n")
    }
    # Gradients, FIM at the new parms
    likobjNew <- likfun_GDFIM(parmsNew)
    if(anyNA(likobjNew)){
      warning(paste("NA found in likobj at iteration", i, "early return...\n"))
      return(list(covparms = parmsNew, loglik = likobjNew$loglik, 
                  grad = likobjNew$grad,
                  info = likobjNew$info, neval = i))
    }
    # Stopping condition
    objNew <- - likobjNew$loglik
    if((abs(sum((parmsNew - parms) * grad)) < convtolOut))
      break
    parms <- parmsNew
    likobj <- likobjNew
    obj <- objNew
  }
  if(!silent)
  {
    cat(paste0("CQCD iter ", i, ": \n"))
    cat("pars = ",  paste0(round(parms, 4)), "  \n" )
    cat(paste0("obj = ", round(obj, 6), "  \n"))
    cat("\n")
  }
  return(list(covparms = parms, loglik = likobj$loglik, loglikInit = loglikInit,
              grad = likobj$grad, info = likobj$info, neval = i + 1,
              betahat = likobj$betahat))
}

#' Coordinate descent using the 2nd order approximation with mini-batching
#' 
#' @param likfun stochastic likelihood function, returns log-likelihood
#' @param likfun_GDFIM returns stochastic log-likelihood, gradient, and FIM
#' @param n the number of total observations
#' @param batchSz batch size
#' @param parms0 starting values of parameters
#' @param convtolOut convergence tolerance on step dot grad
#' @param convtolIn convergence tolerance on the step of one coordinate descent epoch
#' @param cAmij the Armijo condition constant
#' @param maxIterOut maximum number of 2nd order approximations
#' @param maxIterIn maximum number of epochs in coordinate descent
#' @param lb the lower bounds for all parameters
#' @param arg_check the function for checking parms
#' @param silent TRUE/FALSE for suppressing output
CQCD_stochastic <- function(likfun, likfun_GDFIM, n, batchSz, parms0, 
                 convtolOut = 1e-4, convtolIn = 1e-4, 
                 cAmij = 1e-4, maxIterOut = 20, maxIterIn = 40, 
                 lb = rep(0, length(parms0)),
                 arg_check = function(...){T}, silent = FALSE)
{
  if(any(parms0 < lb))
    stop(paste("Initial values for parameters ", which(parms0 < lb), 
             "are smaller than their lower bounds. Stopping...  \n"))
  if(!arg_check(parms0))
    stop(paste("Argument check for parms0 failed. Stopping... \n"))

  parms <- parms0
  batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
  likobj <- likfun_GDFIM(parms, batchIdx)
  
  s <- 0
  scaler <- 1
  iterCounter <- 0
  parmsM1 <- NA
  parmsM2 <- NA
  
  for(i in 1 : maxIterOut)
  {
    obj <- - likobj$loglik
    grad <- - likobj$grad * scaler
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("CQCD iter ", i - 1, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("\n")
    }
    # coordinate descent
    b <- grad - as.vector(H %*% parms)
    DSCObj <- dsc_lb(H, b, parms, silent, convtolIn, maxIterIn, lb)
    parmsNew <- bktk_Armijo(parms, obj, grad, DSCObj$parms, cAmij, 
                            function(x){- likfun(x, batchIdx)$loglik}, 
                            lb, arg_check)
    # if coord descent meets Armijo condition
    if((abs(sum((parmsNew - parms) * grad)) >= convtolOut))
    {
      if(!silent)
        cat("Lower bounded coordinate descent with FIM is used \n")
    }
    else # Use lower bounded gradient descent
    {
      parmsNew <- gd_Armijo(parms, obj, grad, cAmij, 
                            function(x){- likfun(x, batchIdx)$loglik}, 
                            lb, arg_check)
      if(!silent)
        cat("Lower bounded gradient descent is used \n")
    }
    # Stopping condition
    if((abs(sum((parmsNew - parms) * grad)) < convtolOut))
      break
    # Gradients, FIM at the new parms
    batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
    parmsM2 <- parmsM1
    parmsM1 <- parms
    parms <- parmsNew
    iterCounter <- iterCounter + 1
    likobj <- likfun_GDFIM(parms, batchIdx)
    if(i > 1)
      s <- s + sum((parms - parmsM1) * (parmsM1 - parmsM2))
    if(s < 0 && iterCounter > 2){
      s <- 0
      scaler <- scaler / 2
      iterCounter <- 0
      # if(!silent)
      #   cat("scaler is halved\n")
    }
  }
  if(!silent)
  {
    cat(paste0("CQCD iter ", i, ": \n"))
    cat("pars = ",  paste0(round(parms, 4)), "  \n" )
    cat(paste0("obj = ", round(obj, 6), "  \n"))
    cat("\n")
  }
  return(list(covparms = parms, loglik = likobj$loglik,
              grad = likobj$grad, info = likobj$info, neval = i + 1,
              betahat = likobj$betahat))
}

#' Coordinate descent using the 2nd order approximation with mini-batching and
#'   penalty tapering
#' 
#' @param likfun stochastic likelihood function, returns log-likelihood
#' @param likfun_GDFIM returns stochastic log-likelihood, gradient, and FIM
#' @param n the number of total observations
#' @param batchSz batch size
#' @param parms0 starting values of parameters
#' @param convtolOut convergence tolerance on step dot grad
#' @param convtolIn convergence tolerance on the step of one coordinate descent epoch
#' @param cAmij the Armijo condition constant
#' @param maxIterOut maximum number of 2nd order approximations
#' @param maxIterIn maximum number of epochs in coordinate descent
#' @param lb the lower bounds for all parameters
#' @param arg_check the function for checking parms
#' @param silent TRUE/FALSE for suppressing output
CQCD_stochastic_taper <- function(likfun, likfun_GDFIM, n, batchSz, parms0, 
                                  convtolOut = 1e-4, convtolIn = 1e-4, 
                                  cAmij = 1e-4,
                                  maxIterOut = 20, maxIterIn = 40, 
                                  lb = rep(0, length(parms0)),
                                  arg_check = function(...){T}, silent = FALSE)
{
  if(any(parms0 < lb))
    stop(paste("Initial values for parameters ", which(parms0 < lb), 
               "are smaller than their lower bounds. Stopping...  \n"))
  if(!arg_check(parms0))
    stop(paste("Argument check for parms0 failed. Stopping... \n"))
  
  parmsAvg <- parms0
  parms <- parms0
  batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
  likobj <- likfun_GDFIM(parms, batchIdx, parmsAvg, 1)
  
  s <- 0
  scaler <- 1
  iterCounter <- 0
  parmsM1 <- NA
  parmsM2 <- NA
  
  for(i in 1 : maxIterOut)
  {
    obj <- - likobj$loglik
    grad <- - likobj$grad * scaler
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("CQCD iter ", i - 1, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("\n")
    }
    # coordinate descent
    b <- grad - as.vector(H %*% parms)
    DSCObj <- dsc_lb(H, b, parms, silent, convtolIn, maxIterIn, lb)
    parmsNew <- bktk_Armijo(parms, obj, grad, DSCObj$parms, cAmij, 
                            function(x){- likfun(x, batchIdx, 
                                                 (parmsAvg * i + x) / (i + 1),
                                                 i + 1)$loglik}, 
                            lb, arg_check)
    # if coord descent meets Armijo condition
    if((abs(sum((parmsNew - parms) * grad)) >= convtolOut))
    {
      if(!silent)
        cat("Lower bounded coordinate descent with FIM is used \n")
    }
    else # Use lower bounded gradient descent
    {
      parmsNew <- gd_Armijo(parms, obj, grad, cAmij, 
                            function(x){- likfun(x, batchIdx, 
                                                 (parmsAvg * i + x) / (i + 1),
                                                 i + 1)$loglik}, 
                            lb, arg_check)
      if(!silent)
        cat("Lower bounded gradient descent is used \n")
    }
    # Stopping condition
    if((abs(sum((parmsNew - parms) * grad)) < convtolOut))
      break
    # Gradients, FIM at the new parms
    batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
    parmsM2 <- parmsM1
    parmsM1 <- parms
    parms <- parmsNew
    parmsAvg <- (parmsAvg * i + parms) / (i + 1)
    iterCounter <- iterCounter + 1
    likobj <- likfun_GDFIM(parms, batchIdx, parmsAvg, i + 1)
    if(i > 1)
      s <- s + sum((parms - parmsM1) * (parmsM1 - parmsM2))
    if(s < 0 && iterCounter > 2){
      s <- 0
      scaler <- scaler / 2
      iterCounter <- 0
      # if(!silent)
      #   cat("scaler is halved\n")
    }
  }
  if(!silent)
  {
    cat(paste0("CQCD iter ", i, ": \n"))
    cat("pars = ",  paste0(round(parms, 4)), "  \n" )
    cat(paste0("obj = ", round(obj, 6), "  \n"))
    cat("\n")
  }
  return(list(covparms = parms, loglik = likobj$loglik,
              grad = likobj$grad, info = likobj$info, neval = i + 1,
              betahat = likobj$betahat))
}

#' Coordinate descent using the 2nd order approximation with mini-batching and
#'   gradient averaging
#' 
#' @param likfun stochastic likelihood function, returns log-likelihood
#' @param likfun_GDFIM returns stochastic log-likelihood, gradient, and FIM
#' @param n the number of total observations
#' @param batchSz batch size
#' @param parms0 starting values of parameters
#' @param convtolOut convergence tolerance on step dot grad
#' @param convtolIn convergence tolerance on the step of one coordinate descent epoch
#' @param cAmij the Armijo condition constant
#' @param maxIterOut maximum number of 2nd order approximations
#' @param maxIterIn maximum number of epochs in coordinate descent
#' @param lb the lower bounds for all parameters
#' @param avgNum the average window for grad and info
#' @param arg_check the function for checking parms
#' @param silent TRUE/FALSE for suppressing output
CQCD_stochastic_SAG <- function(likfun, likfun_GDFIM, n, batchSz, parms0, 
                            convtolOut = 1e-4, convtolIn = 1e-4, 
                            cAmij = 1e-4, maxIterOut = 20, maxIterIn = 40, 
                            lb = rep(0, length(parms0)), avgNum = 5,
                            arg_check = function(...){T}, silent = FALSE)
{
  if(any(parms0 < lb))
    stop(paste("Initial values for parameters ", which(parms0 < lb), 
               "are smaller than their lower bounds. Stopping...  \n"))
  if(!arg_check(parms0))
    stop(paste("Argument check for parms0 failed. Stopping... \n"))
  
  parms <- parms0
  likobjLst <- list() 
  idxRm <- 1
  for(i in 1 : avgNum){
    batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
    likobjLst[[i]] <- likfun_GDFIM(parms, batchIdx)
  }
  get_likobj <- function(likobjLst, avgNum){
    likobj <- reduce(likobjLst, function(x ,y){list(loglik = x$loglik + y$loglik,
                                                    grad = x$grad + y$grad,
                                                    info = x$info + y$info)})
    likobj$loglik <- likobj$loglik / avgNum
    likobj$grad <- likobj$grad / avgNum
    likobj$info <- likobj$info / avgNum
    likobj
  }
  likobj <- get_likobj(likobjLst, avgNum)
  
  s <- 0
  scaler <- 1
  iterCounter <- 0
  parmsM1 <- NA
  parmsM2 <- NA
  
  for(i in 1 : maxIterOut)
  {
    obj <- - likobj$loglik
    grad <- - likobj$grad * scaler
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("CQCD iter ", i - 1, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("\n")
    }
    # coordinate descent
    b <- grad - as.vector(H %*% parms)
    DSCObj <- dsc_lb(H, b, parms, silent, convtolIn, maxIterIn, lb)
    parmsNew <- bktk_Armijo(parms, obj, grad, DSCObj$parms, cAmij, 
                            function(x){- likfun(x, batchIdx)$loglik}, 
                            lb, arg_check)
    # if coord descent meets Armijo condition
    if((abs(sum((parmsNew - parms) * grad)) >= convtolOut))
    {
      if(!silent)
        cat("Lower bounded coordinate descent with FIM is used \n")
    }
    else # Use lower bounded gradient descent
    {
      parmsNew <- gd_Armijo(parms, obj, grad, cAmij, 
                            function(x){- likfun(x, batchIdx)$loglik}, 
                            lb, arg_check)
      if(!silent)
        cat("Lower bounded gradient descent is used \n")
    }
    # Stopping condition
    if((abs(sum((parmsNew - parms) * grad)) < convtolOut))
      break
    # Gradients, FIM at the new parms
    batchIdx <- sample(x = 1 : n, size = batchSz, replace = F)
    parmsM2 <- parmsM1
    parmsM1 <- parms
    parms <- parmsNew
    iterCounter <- iterCounter + 1

    likobjLst[[idxRm]] <- likfun_GDFIM(parms, batchIdx)
    idxRm <- (idxRm) %% avgNum + 1
    likobj <- get_likobj(likobjLst, avgNum)
    
    if(i > 1)
      s <- s + sum((parms - parmsM1) * (parmsM1 - parmsM2))
    if(s < 0 && iterCounter > 2){
      s <- 0
      scaler <- scaler / 2
      iterCounter <- 0
      # if(!silent)
      #   cat("scaler is halved\n")
    }
  }
  if(!silent)
  {
    cat(paste0("CQCD iter ", i, ": \n"))
    cat("pars = ",  paste0(round(parms, 4)), "  \n" )
    cat(paste0("obj = ", round(obj, 6), "  \n"))
    cat("\n")
  }
  return(list(covparms = parms, loglik = likobj$loglik,
              grad = likobj$grad, info = likobj$info, neval = i + 1,
              betahat = likobj$betahat))
}

#' Coordinate descent for a quadratic function in a lower-bounded domain
#' 
#' @param A 2nd order coefficient matrix (\frac{1}{2} x^\top A x)
#' @param b 1st order coefficient vector
#' @param start_parms starting values of parameters
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of epochs in coordinate descent
#' @param lb the lower bound for the parameters, same length as start_parms
#' 
#' Return a list of two
#' 
#' @return code 0 if convtol is reached, 1 if max number of epochs reached, 2 parms become invalid
#' @return parms new parameter values
dsc_lb <- function(A, b, start_parms, silent, convtol, max_iter, lb)
{
  if(any(diag(A) < 0))
    # stop("A is not non-negative definite. Illegal input for dsc_lb \n")
    return(list(code = -1, parms = rep(0, length(start_parms))))
  if(any(diag(A) == 0 & b < 0))
    # stop("Quadratic will reach Inf. Illegal input for dsc_lb \n")
    return(list(code = -1, parms = rep(0, length(start_parms))))
  parms_new <- start_parms
  optIdx <- setdiff(c(1 : length(start_parms)), which(diag(A) == 0 & b == 0))
  for(k in 1 : max_iter)
  {
    parms_new_cp <- parms_new
    for(j in optIdx)
    {
      if( A[j, j] > 0)
      {
        chg <- - sum(A[j, -j] * parms_new[-j])
        parms_new[j] <- max((- b[j] + chg) / A[j, j], lb[j])
      }else
      {
        parms_new[j] <- lb[j]
      }
      if(is.na(parms_new[j]) || is.infinite(parms_new[j]))
        stop(paste("NA or Inf found in dsc_lb. Early return ...\n"))
    }
    if(sqrt(sum((parms_new - parms_new_cp)^2)) < convtol)
      return(list(code = 0, parms = parms_new))
  }
  return(list(code = 0, parms = parms_new))
}

#' Decide if Armijo condition is satisfied:
#' @param parmsNew the new parameters
#' @param objNew objective function value at parmsNew
#' @param parms the old parameters
#' @param obj objective function value at parms
#' @param grad gradient of the objective function at parms
#' @param c the Armijo rule parameter
#' 
#' @return alpha alpha * d would be the step to take, 
#'   if no proper step size can be found, return -1
is_Armijo <- function(parmsNew, objNew, parms, obj, grad, c)
{
  objNew < obj + c * sum(grad * (parmsNew - parms))
}

#' Optimize with lower bounded gradient descent using Armijo rule
#' @param parms initial parameters
#' @param obj objective function value at parms
#' @param grad gradient at parms0
#' @param c the Armijo rule parameter
#' @param obj_func the objective function that returns a value
#' @param lb the lower bounds of the parameters
#' @param arg_check a function for checking the validity of the input to obj_func
gd_Armijo <- function(parms, obj, grad, c, obj_func, lb, arg_check)
{
  idxRel <- 2 : (length(grad) - 1)
  alphaVec <- parms[idxRel] / grad[idxRel]  
  alphaVec <- alphaVec[alphaVec > 0]
  alphaVec <- alphaVec[is.finite(alphaVec)]
  if(length(alphaVec) == 0)
    alpha <- 1
  else
    alpha <- min(alphaVec)
  for(i in 1 : 7)
  {
    parmsNew <- parms
    parmsNew[idxRel] <- parmsNew[idxRel] - grad[idxRel] * alpha
    parmsNew[parmsNew < lb] <- lb[parmsNew < lb]
    if(arg_check(parmsNew))
    {
      objNew <- obj_func(parmsNew)    
      if(is_Armijo(parmsNew, objNew, parms, obj, grad, c))
        return(parmsNew)
    }
    alpha <- alpha / 4
  }
  warning(paste("gd_Armijo cannot find a new parameter vector that satisfies",
              "the Armijo condition \n"))
  return(parms)
}

#' Backtracking using Armijo rule
#' @param parms old parameters
#' @param obj objective function value at parms
#' @param grad gradient at parms
#' @param parmsNew new parameters
#' @param c the Armijo rule parameter
#' @param obj_func the objective function that returns a value
#' @param lb the lower bounds of the parameters
#' @param arg_check a function for checking the validity of the input to obj_func
bktk_Armijo <- function(parms, obj, grad, parmsNew, c, obj_func, lb, arg_check)
{
  idxRel <- 2 : (length(grad) - 1)
  alpha <- 1 
  step <- parmsNew - parms
  while(alpha > 4^(-7))
  {
    parmsTmp <- parms + step * alpha
    parmsTmp[parmsTmp < lb] <- lb[parmsTmp < lb]
    if(arg_check(parmsTmp))
    {
      objNew <- obj_func(parmsTmp)    
      if(is_Armijo(parmsTmp, objNew, parms, obj, grad, c))
        return(parmsTmp)
    }
    alpha <- alpha / 4
  }
  return(parms)
}

#' Wrapper for different flavors of CQCD
#'
#' @param theta starting values of parameters 
#' @param y obs vec
#' @param locs covariate matrix
#' @param NNarray NNarray
#' @param lb the lower bounds for all parameters
#' @param lambda penalty parm
#' @param pen_fun penalty func
#' @param dpen_fun dpenalty func
#' @param ddpen_fun ddpenalty func
#' @param convCQCD convergence level of CQCD
#' @param convCCD convergence level of CCD
#' @param cAmij Armijo constant
#' @param maxIterCQCD max # iterations of CQCD
#' @param maxIterCCD max # iterations of CCD
#' @param covFn covariance function name
#' @param miniCQCD mini-batch size for CQCD
#' @param mini bool for using mini-batching?
#' @param taper bool for using penalty tapering?
#' @param silent bool for silent execution?
CQCD_wrap <- function(theta, y, locs, NNarray, lb, lambda, 
                      pen_fun, dpen_fun, ddpen_fun, 
                      convCQCD = 1e-4, convCCD = 1e-4, cAmij = 1e-4, 
                      maxIterCQCD = 500, maxIterCCD = 40,
                      covFn = "matern25_scaledim_sqrelevance",
                      miniCQCD = 128, mini = T, taper = F, silent = F)
{
  n <- nrow(locs)
  if(mini){
    if(taper){
      objfun <- function(theta, batchIdx, thetaAvg, iter){
        likObj <- vecchia_meanzero_loglik(theta, covFn, y, locs, 
                                          NNarray[batchIdx, , drop = F])
        likObj$loglik <- likObj$loglik - pen_fun(thetaAvg, lambda, iter)
        likObj
      }
      objfun_gdfm <- function(theta, batchIdx, thetaAvg, iter){
        likObj <- vecchia_meanzero_loglik_grad_info(theta, covFn, y, locs,
                                                    NNarray[batchIdx, , 
                                                            drop = F])
        likObj$loglik <- likObj$loglik - pen_fun(thetaAvg, lambda, iter)
        likObj$grad <- likObj$grad - dpen_fun(thetaAvg, lambda, iter)
        likObj$info <- likObj$info + ddpen_fun(thetaAvg, lambda, iter)
        likObj
      }
      CQCD_stochastic_taper(objfun, objfun_gdfm, n, miniCQCD, 
                            theta, lb = lb, convtolOut = convCQCD, 
                            convtolIn = convCCD, maxIterOut = maxIterCQCD,
                            maxIterIn = maxIterCCD, cAmij = cAmij, 
                            silent = silent)
    }else{
      objfun <- function(theta, batchIdx){
        likObj <- vecchia_meanzero_loglik(theta, covFn, y, locs, 
                                          NNarray[batchIdx, , drop = F])
        likObj$loglik <- likObj$loglik - pen_fun(theta, lambda)
        likObj
      }
      objfun_gdfm <- function(theta, batchIdx){
        likObj <- vecchia_meanzero_loglik_grad_info(theta, covFn, y, locs,
                                                    NNarray[batchIdx, , 
                                                            drop = F])
        likObj$loglik <- likObj$loglik - pen_fun(theta, lambda)
        likObj$grad <- likObj$grad - dpen_fun(theta, lambda)
        likObj$info <- likObj$info + ddpen_fun(theta, lambda)
        likObj
      }
      CQCD_stochastic(objfun, objfun_gdfm, n, miniCQCD, 
                      theta, lb = lb, convtolOut = convCQCD, 
                      convtolIn = convCCD, maxIterOut = maxIterCQCD,
                      maxIterIn = maxIterCCD, cAmij = cAmij, 
                      silent = silent)
    }
  }else{
    stop("CQCD_wrap with mini = F is not implemented yet\n")
  }
}
