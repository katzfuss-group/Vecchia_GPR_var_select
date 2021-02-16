library(GpGp)

#' Coordinate descent using the 2nd order approximation of the log-likelihood 
#' 
#' @param likfun likelihood function, returns log-likelihood
#' @param likfunGDFIM returns log-likelihood, gradient, and FIM
#' @param locs the n-by-d location matrix 
#' @param p the number of columns in the design matrix X
#' @param start_parms starting values of parameters
#' @param lambda L1 penalty parameter
#' @param epsl step size when coordinate descent does not reduce the obj func
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on step dot grad
#' @param convtol2 convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of 2nd order approximations
#' @param max_iter2 maximum number of epochs in coordinate descent
#' @param lb_nonrel_parms the lower bounds for non-relevance parameters
quad_cdsc_L1 <- function(likfun, likfunGDFIM, locs, p, start_parms, lambda, 
                         epsl, silent = FALSE, convtol = 1e-4, 
                         convtol2 = 1e-4, max_iter = 40, max_iter2 = 40, 
                         lb_nonrel_parms = rep(0, length(start_parms) - ncol(locs)))
{
  if(lambda < 0 || epsl < 0)
    stop("lambda and epsl should both be greater than zero\n")
  lb_parms <- c(lb_nonrel_parms[1], rep(0, ncol(locs)), lb_nonrel_parms[-1])
  if(any(start_parms < lb_parms))
    stop("some coefficient in start_parms is smaller than its lower bound")
  parms <- start_parms
  nloc <- ncol(locs)
  idxPosiLocs <- parms[2 : (1 + nloc)] > 0
  idxPosiParm <- c(T, idxPosiLocs, rep(T, length(parms) - 1 - nloc))
  parmsPosi <- parms[idxPosiParm]
  # check if the relevance parameters are all too close to zero
  if(sum(idxPosiLocs) == 0)
  {
    cat("quad_cdsc_L1 reached zero for all relevance parameters\n")
    return(list(covparms = parms, obj = NA, betahat = rep(0, p), neval = 0))
  }
  likobj <- likfunGDFIM(parmsPosi, locs[, idxPosiLocs])
  for(i in 1 : max_iter)
  {
    # check for Inf, NA, or NaN
    if( !GpGp::test_likelihood_object(likobj) ){
      stop("inf or na or nan in likobj\n")
    }
      
    obj <- -likobj$loglik + lambda * sum(parmsPosi[2 : (1 + sum(idxPosiLocs))])
    grad <- -likobj$grad 
    grad[2 : (1 + sum(idxPosiLocs))] <- grad[2 : (1 + sum(idxPosiLocs))] + lambda
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("Iter ", i, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("\n")
    }
    # coordinate descent
    b <- grad - as.vector(H %*% parmsPosi)
    coord_des_obj <- cdsc_quad_posi(H, b, parmsPosi, silent, 
                                    convtol2, max_iter2, 1e6, lb_parms[idxPosiParm])
    # check if the relevance parameters are all zero
    if(sum(coord_des_obj$parms[2 : (1 + sum(idxPosiLocs))]) == 0)
      coord_des_obj$code <- 2
    
    # check if obj func decreases
    if(coord_des_obj$code < 2) # parms_new is valid
    {
      stepSz <- step_Armijo(parmsPosi, obj, grad, coord_des_obj$parms - parmsPosi, 1e-4, 
                            function(x){- likfun(x, locs[, idxPosiLocs])$loglik + 
                                lambda * sum(x[2 : (1 + sum(idxPosiLocs))])})
      if(stepSz < 0)
        grad_des <- T
      else
      {
        parmsNew <- parms
        parmsNew[idxPosiParm] <- parmsNew[idxPosiParm] + stepSz * (coord_des_obj$parms - parmsPosi)
        grad_des <- F
      }
    }
    else
      grad_des <- T
    
    if(grad_des)
    {
      parmsNewPosiTmp <- step_grad_Armijo(parms[idxPosiParm], c(2 : (sum(idxPosiParm) - 1)), 
                                          obj, grad, lb_parms[idxPosiParm], 1e-3, 
                                       function(x){- likfun(x, locs[, idxPosiLocs])$loglik + 
                                           lambda * sum(x[2 : (1 + sum(idxPosiLocs))])},
                                       function(x){sum(x[2 : (1 + sum(idxPosiLocs))]) > 1e-3})
      if(all(parmsNewPosiTmp == -1))
      {
        cat("gradient descent cannot find proper parameter values\n")
        return(list(covparms = lb_parms, obj = NA, betahat = rep(0, p), neval = i))
      }
      
      parmsNew <- parms
      parmsNew[idxPosiParm] <- parmsNewPosiTmp
      cat("Gradient descent is used\n")
    }
    
    if(!silent)
      cat("\n")
    
    idxPosiLocs <- parmsNew[2 : (1 + nloc)] > 0
    idxPosiParmNew <- c(T, idxPosiLocs, rep(T, length(parmsNew) - 1 - nloc))
    parmsNewPosi <- parmsNew[idxPosiParmNew]
    likobjNew <- likfunGDFIM(parmsNewPosi, locs[, idxPosiLocs])
    objNew <- - likobjNew$loglik + lambda * sum(parmsNewPosi[2 : (1 + sum(idxPosiLocs))])
    if((objNew > obj) || 
       (abs(sum((parmsNew - parms)[idxPosiParm] * (grad))) < convtol))
      break
    idxPosiParm <- idxPosiParmNew
    parms <- parmsNew
    parmsPosi <- parmsNewPosi
    likobj <- likobjNew
    obj <- objNew
  }
  return(list(covparms = parms, obj = obj, betahat = likobj$betahat, neval = i + 1))
}

#' Coordinate descent for a quadratic function in the positive domain 
#' 
#' @param A 2nd order coefficient matrix (\frac{1}{2} x^\top A x)
#' @param b 1st order coefficient vector
#' @param start_parms starting values of parameters
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of epochs in coordinate descent
#' @param max_parm maximum parameter value
#' @param lb the lower bound for the parameters, same length as start_parms
#' 
#' Return a list of two
#' 
#' @return code 0 if convtol is reached, 1 if max number of epochs reached, 2 parms become invalid
#' @return parms new parameter values
cdsc_quad_posi <- function(A, b, start_parms, silent, convtol, max_iter, max_parm, lb)
{
  parms_new <- start_parms
  for(k in 1 : max_iter)
  {
    parms_new_cp <- parms_new
    for(j in 1 : length(parms_new))
    {
      chg <- - sum(A[j, -j] * parms_new[-j])
      parms_new[j] <- max((- b[j] + chg) / A[j, j], lb[j])
      if(parms_new[j] > max_parm || is.na(parms_new[j]) || is.infinite(parms_new[j]))
        return(list(code = 2, parms = parms_new))
    }
    if(sqrt(sum((parms_new - parms_new_cp)^2)) < convtol)
      return(list(code = 0, parms = parms_new))
  }
  return(list(code = 1, parms = parms_new))
}

#' Find step size using Armijo rule:
#'   the initial step size is assumed one
#'   half the step size each iteration until Armijo rule is satisfied
#' @param parms0 initial parameters
#' @param v0 objective function value at parms0
#' @param d0 gradient of the objective function at parms0
#' @param d a descent direction (norm 1 not required)
#' @param c the Armijo rule parameter
#' @param obj_func the objective function that returns a value
#' 
#' @return alpha alpha * d would be the step to take, 
#'   if no proper step size can be found, return -1
step_Armijo <- function(parms0, v0, d0, d, c, obj_func)
{
  alpha <- 1
  step <- alpha * d
  val <- obj_func(parms0 + step)
  while(val > v0 + c * sum(d0 * step))
  {
    alpha <- 0.5 * alpha
    step <- alpha * d
    if(max(abs(step)) < 1e-5)
      return(-1)
    val <- obj_func(parms0 + step)
  }
  return(alpha)
}

#' Find step size based on gradient using Armijo rule:
#' @param parms0 initial parameters
#' @param idxRel indices for the relevance parameters
#' @param v0 objective function value at parms0
#' @param d0 gradient of the objective function at parms0
#' @param lb the lower bounds of the parameters
#' @param c the Armijo rule parameter
#' @param obj_func the objective function that returns a value
#' @param arg_check a function for checking the validity of the input to obj_func
#' 
#' @return parms if valid parms is found, otherwise -1
step_grad_Armijo <- function(parms0, idxRel, v0, d0, lb, c, obj_func, arg_check)
{
  alphaVec <- sort(((parms0 - lb) / d0)[idxRel], decreasing = T)
  alphaVec <- alphaVec[alphaVec > 0]
  for(alpha in alphaVec)
  {
    parms <- parms0
    parms[idxRel] <- parms0[idxRel] - alpha * d0[idxRel]
    parms[parms < lb] <- lb[parms < lb]
    if(arg_check(parms) == F)
      next
    step <- parms - parms0
    val <- obj_func(parms)
    if(val < v0 + c * sum(d0 * step))
      return(parms)
  }
  return(-1)
}