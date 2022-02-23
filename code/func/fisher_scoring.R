source("misc.R")

fisher_scoring_meanzero <- function (likfun, start_parms, link, silent = FALSE, 
                                     convtol = 1e-04, max_iter = 40) 
{
  wolfe_check <- function(likobj0, likobj1, logparms, step, 
                          both) {
    c1 <- 1e-04
    c2 <- 0.9
    tol <- 0.1
    ll0 <- likobj0$loglik
    gr0 <- likobj0$grad
    ll1 <- likobj1$loglik
    gr1 <- likobj1$grad
    if (!both) {
      satfd <- ll1 <= ll0 + c1 * crossprod(step, gr0) + 
        tol
    }
    else {
      satfd <- ll1 <= ll0 + c1 * crossprod(step, gr0) + 
        tol && -crossprod(step, gr1) <= -c2 * crossprod(step, 
                                                        gr0) + tol
    }
    return(satfd)
  }
  logparms <- start_parms
  likobj <- likfun(logparms)
  if (!test_likelihood_object(likobj)) {
    logparms <- 0.1 * logparms
    likobj <- likfun(logparms)
  }
  loglik <- likobj$loglik
  grad <- likobj$grad
  info <- likobj$info
  diag(info) <- diag(info) + 0.1 * min(diag(info))
  if (!silent) {
    cat(paste0("Iter ", 0, ": \n"))
    cat("pars = ", paste0(round(link(logparms), 4)), "  \n")
    cat(paste0("loglik = ", round(-loglik, 6), "  \n"))
    cat("grad = ")
    cat(as.character(round(-grad, 3)))
    cat("\n\n")
  }
  for (j in 1:max_iter) {
    likobj0 <- likobj
    tol <- 1e-10
    if (condition_number(info) > 1/tol) {
      if (!silent) 
        cat("Cond # of info matrix > 1/tol \n")
      diag(info) <- diag(info) + tol * max(diag(info))
    }
    step <- -solve(info, grad)
    if (mean(step^2) > 1) {
      if (!silent) 
        cat("##\n")
      step <- step/sqrt(mean(step^2))
    }
    newlogparms <- logparms + step
    likobj <- likfun(newlogparms)
    cnt <- 1
    while (!test_likelihood_object(likobj)) {
      if (!silent) 
        cat("inf or na or nan in likobj\n")
      step <- 0.5 * step
      newlogparms <- logparms + step
      likobj <- likfun(newlogparms)
      if (cnt == 10) {
        stop("could not avoid inf, na or nan\n")
      }
    }
    cnt <- 1
    no_decrease <- FALSE
    both <- FALSE
    mult <- 1
    stepgrad <- c(crossprod(step, grad))
    logparms <- logparms + step
    loglik <- likobj$loglik
    grad <- likobj$grad
    info <- likobj$info
    if (!silent) {
      cat(paste0("Iter ", j, ": \n"))
      cat("pars = ", paste0(round(link(logparms), 4)), 
          "  \n")
      cat(paste0("loglik = ", round(-loglik, 6), "  \n"))
      cat("grad = ")
      cat(as.character(round(grad, 4)), "\n")
      cat("step dot grad = ", stepgrad, "\n")
      cat("\n")
    }
    if (abs(stepgrad) < convtol || no_decrease) {
      break
    }
  }
  # betahatinfo <- likobj
  # betahat <- as.vector(betahatinfo$betahat)
  # betacov <- solve(betahatinfo$betainfo)
  # sebeta <- sqrt(diag(betacov))
  # tbeta <- betahat/sebeta
  ret <- list(covparms = link(logparms), logparms = logparms, 
              # betahat = betahat, sebeta = sebeta, betacov = betacov, 
              # tbeta = tbeta, 
              loglik = loglik, no_decrease = no_decrease, 
              grad = likobj$grad, 
              info = likobj$info, 
              conv = (abs(stepgrad) < 
                        convtol || no_decrease), neval = j + 1)
  return(ret)
}

#' Train model indicated by locsIdx using Fisher scoring 
#' 
#' @param locsIdx covariate idx
#' @param m NN size
#' @param theta starting values of parameters 
#' @param locs in-sample covariate matrix
#' @param y in-sample obs vec
#' @param covFn covariance function name
#' @param maxIter max # iterations of Fisher scoring
#' @param conv convergence level of Fisher scoring
#' @param link_func link func, e.g., log
#' @param resp_func response func, e.g., exp
#' @param dresp_func dresponse func, e.g., exp
#' @param silent bool for silent execution?
fisher_scoring_meanzero_wrap <- function(locsIdx, m, theta, locs, y,
                                         covFn, maxIter = 100, 
                                         conv = 1e-4, link_func = log, 
                                         resp_func = exp, dresp_func = exp, 
                                         silent = T)
{
  d <- ncol(locs)
  crtIter <- 1
  ttlIter <- 0
  while(TRUE){
    # MM and NN
    # use dummy lb (-Inf) for MM_NN
    MMNNObj <- MM_NN(theta, locsIdx, locs, y, rep(-Inf, length(theta)), m)
    # obj func
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
    # Fisher scoring
    FisherObj <- fisher_scoring_meanzero(objfun, link_func(MMNNObj$thetaRel), 
                                         resp_func, silent, conv, crtIter)
    # Update theta
    theta <- rep(0, d + 2)
    theta[c(1, locsIdx + 1, d + 2)] <- resp_func(FisherObj$logparms)
    # Update iter num
    ttlIter <- ttlIter + crtIter
    crtIter <- min(2 * crtIter, maxIter - ttlIter)
    if(ttlIter == maxIter)
      break
  }
  # return
  theta
}






