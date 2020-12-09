library(GpGp)

#' Coordinate descent using the 2nd order approximation of the log-likelihood 
#' 
#' @param likfun likelihood function, returns log-likelihood, gradient, and FIM
#' @param start_parms starting values of parameters
#' @param lambda L1 penalty parameter
#' @param pen_idx the indices of parameters that will be penalized by the L1
#' @param epsl step size when coordinate descent does not reduce the obj func
#' @param silent TRUE/FALSE for suppressing output
#' @param convtol convergence tolerance on the objective function
#' @param convtol2 convergence tolerance on the step of one coordinate descent epoch
#' @param max_iter maximum number of 2nd order approximations
#' @param max_iter2 maximum number of epochs in coordinate descent
quad_cdsc_L1 <- function(likfun, start_parms, lambda, pen_idx, epsl, 
                              silent = FALSE, convtol = 1e-4, convtol2 = 1e-4, 
                              max_iter = 40, max_iter2 = 40)
{
  
  if(lambda < 0 || epsl < 0)
    stop("lambda and epsl should both be greater than zero\n")
  parms <- start_parms
  likobj <- likfun(parms)
  for(i in 1 : max_iter)
  {
    # check for Inf, NA, or NaN
    if( !GpGp::test_likelihood_object(likobj) ){
      stop("inf or na or nan in likobj\n")
    }
    obj <- -likobj$loglik + lambda * sum(abs(parms[pen_idx]))
    grad <- -likobj$grad 
    grad[pen_idx] <- grad[pen_idx] + lambda
    H <- likobj$info
    if(!silent)
    {
      cat(paste0("Iter ", i, ": \n"))
      cat("pars = ",  paste0(round(parms, 4)), "  \n" )
      cat(paste0("obj = ", round(obj, 6), "  \n"))
      cat("grad = ")
      cat(as.character(round(grad, 3)))
      cat("\n")
    }
    
    # coordinate descent
    b <- grad - as.vector(H %*% parms)
    coord_des_obj <- cdsc_quad_posi(H, b, parms, silent, 
                                            convtol2, max_iter2, 1e6)
    # check if obj func decreases
    if(coord_des_obj$code < 2) # parms_new is valid
    {
      stepSz <- step_Armijo(parms, obj, grad, coord_des_obj$parms - parms, 1e-4, 
                                 function(x){- likfun(x)$loglik + lambda * sum(x)})
      if(stepSz < 0)
        grad_des <- T
      else
      {
        parmsNew <- parms + stepSz * (coord_des_obj$parms - parms)
        grad_des <- F
      }
    }
    else
      grad_des <- T
    
    if(grad_des)
    {
      parmsNew <- parms - grad * epsl
      parmsNew[parmsNew < 0] <- 0
      if(!silent)
      {
        cat("Gradient descent is used\n")
      }
    }
    
    if(!silent)
      cat("\n")
    
    likobjNew <- likfun(parmsNew)
    objNew <- - likobjNew$loglik + lambda * sum(parmsNew)
    if(objNew > obj - convtol)
      break
    parms <- parmsNew
    likobj <- likobjNew
  }
  return(list(covparms = parms))
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
#' 
#' Return a list of two
#' 
#' @return code 0 if convtol is reached, 1 if max number of epochs reached, 2 parms become invalid
#' @return parms new parameter values
cdsc_quad_posi <- function(A, b, start_parms, silent, convtol, max_iter, max_parm)
{
  parms_new <- start_parms
  for(k in 1 : max_iter)
  {
    parms_new_cp <- parms_new
    for(j in 1 : length(parms_new))
    {
      chg <- - sum(A[j, -j] * parms_new[-j])
      parms_new[j] <- max((- b[j] + chg) / A[j, j], 0)
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