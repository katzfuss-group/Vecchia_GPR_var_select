# Experiment 2
n <- 1e4
d <- 1e2
r <- c(10, 5, 2, 1, 0.5, rep(0, d - 5))
sigmasq <- 1.0 # variance
tausq <- 0.05^2 # nugget

set.seed(123)
locs <- lhs::randomLHS(n, d)
locs <- locs * outer(rep(sqrt(n), n), 1 / sqrt(colSums(locs^2)))
covM <- GpGp::matern25_scaledim_relevance(c(sigmasq, r, tausq), locs)
cholM <- t(chol(covM))
y <- as.vector(cholM %*% rnorm(n))
X <- matrix(1, n, 1)

rInit <- rep(1, d)
sigmasqInit <- 0.25
tausqInit <- 0.01^2
lambdaVec <- c(0, 1, 10, 100)
lb_nonrel_parms <- c(0.01^2, 0.01^2)
cvInnerIter <- 20
cvOuterIter <- 2

# quad_cdsc_L1 method
theta <- c(sigmasqInit, rInit, tausqInit)
crtIter <- 1
maxIter <- 100
m <- 30
sink("quad_cdsc_L1.out")
# cross-validation
for(i in 1 : cvOuterIter)
{
  # CV inputs
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  est_func <- function(lambda, locs, X, y, NNarray)
  {
    objfun <- function(theta, locs){
      cat("v\n")
      GpGp::vecchia_profbeta_loglik(theta, "matern25_scaledim_relevance",
                                    y, X, locs, NNarray)
    }
    objfun_gdfm <- function(theta, locs){
      cat("g\n")
      GpGp::vecchia_profbeta_loglik_grad_info(theta, "matern25_scaledim_relevance",
                                              y, X, locs, NNarray)
    }
    # Notice that theta will be inherited from the outerloop
    quad_cdsc_L1(objfun, objfun_gdfm, locs, 1, theta, lambda, 1e-3, silent = T, 
                 max_iter = cvInnerIter, max_iter2 = 40, lb_nonrel_parms = lb_nonrel_parms)
  }
  pred_func <- function(rslt, locs_pred, X_pred, locs_obs, X_obs, y_obs)
  {
    idxPosiSub <- rslt$covparms[2 : (length(rslt$covparms) - 1)] > 0
    idxPosi <- c(T, idxPosiSub, T)
    rslt$covparms <- rslt$covparms[idxPosi]
    if(sum(idxPosiSub) == 0)
      return(c(X_pred %*% rslt$betahat))
    GpGp::predictions(fit = rslt, locs_pred = locs_pred[, idxPosiSub], X_pred = X_pred, 
                      y_obs = y_obs, locs_obs = locs_obs[, idxPosiSub], X_obs = X_obs, 
                      covfun_name = "matern25_scaledim_relevance")
  }
  # Compute CV loss
  loss <- rep(NA, length(lambdaVec))
  idxRnd <- sample(c(1 : n), n, F)
  for(j in 1 : length(lambdaVec))
  {
    loss[j] <- cross_valid(est_func, pred_func, crit_MSE, locsOdr, XOdr, yOdr, 
                           m, lambdaVec[j], 5, idxRnd)
    cat(lambdaVec[j], ": loss ", loss[j], "\n")
  }
  lambda <- lambdaVec[which.min(loss)]
  cat("CV", i, ": best lambda = ", lambda, "\n")
  # Fit with all data
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  theta <- est_func(lambda, locsOdr, XOdr, yOdr, NNarray)$covparms
  cat("CV", i, ": parms = ", theta, "\n")
}
# final fitting
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  # Define functions for parameter estimation in the outer loop
  objfun <- function(theta, locsOdr){
    cat("v\n")
    GpGp::vecchia_profbeta_loglik(theta, "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
  }
  objfun_gdfm <- function(theta, locsOdr){
    cat("g\n")
    GpGp::vecchia_profbeta_loglik_grad_info(theta, "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
  }
  
  theta <- quad_cdsc_L1(objfun, objfun_gdfm, locsOdr, 1, theta, lambda, 1e-3, silent = T, 
               max_iter = min(crtIter, maxIter - crtIter + 1), max_iter2 = 40, 
               lb_nonrel_parms = lb_nonrel_parms)$covparms
  cat("quad_cdsc_L1 fit iter", crtIter, ": estimated parms = ", theta, "\n")
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# Fisher scoring
tausqInitFS <- tausqInit + 0.01^2
theta <- c(sigmasqInit, rInit, tausqInitFS)
crtIter <- 1
maxIter <- 100
m <- 30

linkfuns <- GpGp::get_linkfun("matern25_scaledim")
link <- linkfuns$link
dlink <- linkfuns$dlink
invlink <- linkfuns$invlink
thetaTrans <- invlink(theta)

pen <- function(theta){lambda * sum(theta[2 : (d + 1)])}
dpen <- function(theta){lambda * c(0, rep(1, d), 
                                   rep(0, length(theta) - d - 1))}
ddpen <- function(theta){matrix(0, length(theta), length(theta))}
sink("FS_relevance.out")
# cross-validation
for(i in 1 : cvOuterIter)
{
  # CV inputs
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  est_func <- function(lambda, locs, X, y, NNarray)
  {
    objfun <- function(thetaTrans){
      cat("g")
      likobj <- 
        GpGp::vecchia_profbeta_loglik_grad_info(link(thetaTrans), 
                                                "matern25_scaledim_relevance",
                                                y, X, locs, NNarray)
      likobj$loglik <- -likobj$loglik + pen(link(thetaTrans))
      likobj$grad <- -c(likobj$grad)*dlink(thetaTrans) +
        dpen(link(thetaTrans))*dlink(thetaTrans)
      likobj$info <- likobj$info*outer(dlink(thetaTrans),dlink(thetaTrans)) +
        ddpen(link(thetaTrans))*outer(dlink(thetaTrans),dlink(thetaTrans))
      return(likobj)
    }
    # Notice that theta will be inherited from the outerloop
    fisher_scoring(objfun, thetaTrans, link, T, 1e-3, cvInnerIter)
  }
  pred_func <- function(rslt, locs_pred, X_pred, locs_obs, X_obs, y_obs)
  {
    rslt$covparms <- link(rslt$logparms)
    idxPosiSub <- rslt$covparms[2 : (length(rslt$covparms) - 1)] > 0
    idxPosi <- c(T, idxPosiSub, T)
    rslt$covparms <- rslt$covparms[idxPosi]
    if(sum(idxPosiSub) == 0)
      return(c(X_pred %*% rslt$betahat))
    GpGp::predictions(fit = rslt, locs_pred = locs_pred[, idxPosiSub], X_pred = X_pred, 
                      y_obs = y_obs, locs_obs = locs_obs[, idxPosiSub], X_obs = X_obs, 
                      covfun_name = "matern25_scaledim_relevance")
  }
  # Compute CV loss
  loss <- rep(NA, length(lambdaVec))
  idxRnd <- sample(c(1 : n), n, F)
  for(j in 1 : length(lambdaVec))
  {
    lambda <- lambdaVec[j]
    loss[j] <- cross_valid(est_func, pred_func, crit_MSE, locsOdr, XOdr, yOdr, 
                           m, lambdaVec[j], 5, idxRnd)
    cat(lambdaVec[j], ": loss ", loss[j], "\n")
  }
  lambda <- lambdaVec[which.min(loss)]
  cat("CV", i, ": best lambda = ", lambda, "\n")
  # Fit with all data
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  thetaTrans <- est_func(lambda, locsOdr, XOdr, yOdr, NNarray)$logparms
  cat("CV", i, ": parms = ", theta, "\n")
}
# final fitting
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(link(thetaTrans[2 : (d + 1)]))
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)

  objfun <- function(thetaTrans){
    cat("g")
    likobj <- 
      GpGp::vecchia_profbeta_loglik_grad_info(link(thetaTrans), 
                                              "matern25_scaledim_relevance",
                                              yOdr, XOdr, locsOdr, 
                                              NNarray)
    likobj$loglik <- -likobj$loglik + pen(link(thetaTrans))
    likobj$grad <- -c(likobj$grad)*dlink(thetaTrans) +
      dpen(link(thetaTrans))*dlink(thetaTrans)
    likobj$info <- likobj$info*outer(dlink(thetaTrans),dlink(thetaTrans)) +
      ddpen(link(thetaTrans))*outer(dlink(thetaTrans),dlink(thetaTrans))
    return(likobj)
  }
  thetaTrans <- fisher_scoring(objfun, thetaTrans, link, T, 1e-3,
                        min(crtIter, maxIter - crtIter + 1))$logparms
  cat("Fisher scoring fit iter", crtIter, ": estimated parms = ", link(thetaTrans), "\n")
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# optim package from R
theta <- c(sigmasqInit, rInit, tausqInit)
crtIter <- 1
maxIter <- 100
m <- 30
# sink("L-BFGS-B.out")
# Cross-validation
for(i in 1 : cvOuterIter)
{
  # CV inputs
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  est_func <- function(lambda, locs, X, y, NNarray)
  {
    objfun <- function(theta){
      cat("v\n")
      likobj <- GpGp::vecchia_profbeta_loglik(theta, "matern25_scaledim_relevance",
                                    y, X, locs, NNarray)
      return(-likobj$loglik + lambda * sum(theta[2 : (d + 1)]))
    }
    objfun_gdfm <- function(theta){
      cat("g\n")
      likobj <- GpGp::vecchia_profbeta_loglik_grad_info(theta, "matern25_scaledim_relevance",
                                                        y, X, locs, NNarray)
      return(-likobj$grad + lambda * c(0, rep(1, d), 
                                       rep(0, length(theta) - d - 1)))
    }
    # Notice that theta will be inherited from the outerloop
    thetaOpt <- optim(theta, objfun, objfun_gdfm, method = "L-BFGS-B",
                      lower = c(lb_nonrel_parms[1], rep(0, d), lb_nonrel_parms[-1]), 
                      control = list(maxit = cvInnerIter,
                                     trace = 0,
                                     pgtol = 1e-3))$par
    # compute betahat separately
    likobj <- GpGp::vecchia_profbeta_loglik(thetaOpt, "matern25_scaledim_relevance",
                                            y, X, locs, NNarray)
    return(list(covparms = thetaOpt, betahat = likobj$betahat))
  }
  pred_func <- function(rslt, locs_pred, X_pred, locs_obs, X_obs, y_obs)
  {
    idxPosiSub <- rslt$covparms[2 : (length(rslt$covparms) - 1)] > 0
    idxPosi <- c(T, idxPosiSub, T)
    rslt$covparms <- rslt$covparms[idxPosi]
    if(sum(idxPosiSub) == 0)
      return(c(X_pred %*% rslt$betahat))
    GpGp::predictions(fit = rslt, locs_pred = locs_pred[, idxPosiSub], X_pred = X_pred, 
                      y_obs = y_obs, locs_obs = locs_obs[, idxPosiSub], X_obs = X_obs, 
                      covfun_name = "matern25_scaledim_relevance")
  }
  # Compute CV loss
  loss <- rep(NA, length(lambdaVec))
  idxRnd <- sample(c(1 : n), n, F)
  for(j in 1 : length(lambdaVec))
  {
    loss[j] <- cross_valid(est_func, pred_func, crit_MSE, locsOdr, XOdr, yOdr, 
                           m, lambdaVec[j], 5, idxRnd)
    cat(lambdaVec[j], ": loss ", loss[j], "\n")
  }
  lambda <- lambdaVec[which.min(loss)]
  cat("CV", i, ": best lambda = ", lambda, "\n")
  # Fit with all data
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  theta <- est_func(lambda, locsOdr, XOdr, yOdr, NNarray)$covparms
  cat("CV", i, ": parms = ", theta, "\n")
}
# Fit model with chosen lambda
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  
  objfun <- function(theta){
    cat("v", "\n")
    likobj <- GpGp::vecchia_profbeta_loglik(theta, 
                                            "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
    return(-likobj$loglik + lambda * sum(theta[2 : (d + 1)]))
  }
  dobjfun <- function(theta){
    cat("g", "\n")
    likobj <- GpGp::vecchia_profbeta_loglik_grad_info(theta, 
                                                      "matern25_scaledim_relevance",
                                                      yOdr, XOdr, locsOdr, 
                                                      NNarray)
    return(-likobj$grad + lambda * c(0, rep(1, d), 
                                     rep(0, length(theta) - d - 1)))
  }
  
  theta <- optim(theta, objfun, dobjfun, method = "L-BFGS-B",
                 lower = c(lb_nonrel_parms[1], rep(0, d), lb_nonrel_parms[-1]), 
                 control = list(maxit = min(crtIter, maxIter - crtIter + 1),
                                trace = 0,
                                pgtol = 1e-3))$par
  cat("Optim fit iter", crtIter, ": estimated parms = ", theta, "\n")
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# Forward selection
lambda <- 0
sink("forward.out")
stop_cond <- function(optObjOld, optObjNew)
{
  BICOld <- log(n) * length(optObjOld$parms) + 2 * optObjOld$obj
  BICNew <- log(n) * length(optObjNew$parms) + 2 * optObjNew$obj
  return(BICNew > BICOld)
}
criteria <- function(varIdx)
{
  tausqInitFW <- tausqInit + 0.01^2
  theta <- c(sigmasqInit, rInit[varIdx], tausqInitFW)
  crtIter <- 1
  maxIter <- 100
  m <- 30
  d <- length(varIdx)
  linkfuns <- get_linkfun("matern25_scaledim")
  link <- linkfuns$link
  dlink <- linkfuns$dlink
  invlink <- linkfuns$invlink
  thetaTrans <- invlink(theta)
  pen <- function(theta){lambda * sum(theta[2 : (d + 1)])}
  dpen <- function(theta){lambda * c(0, rep(1, d), 
                                     rep(0, length(theta) - d - 1))}
  ddpen <- function(theta){matrix(0, length(theta), length(theta))}
  while(maxIter >= crtIter)
  {
    locsScal <- locs[, varIdx, drop = F] %*% diag(x = theta[2 : (length(varIdx) + 1)], 
                                                  nrow = length(varIdx))
    odr <- GpGp::order_maxmin(locsScal)
    yOdr <- y[odr]
    locsOdr <- locs[odr, varIdx, drop = F]
    XOdr <- X[odr, , drop = F]
    NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
    
    objfun1 <- function(thetaTrans){
      cat("g\n")
      likobj <- 
        GpGp::vecchia_profbeta_loglik_grad_info(link(thetaTrans), 
                                                "matern25_scaledim_relevance",
                                                yOdr, XOdr, locsOdr, 
                                                NNarray)
      likobj$loglik <- -likobj$loglik + pen(link(thetaTrans))
      likobj$grad <- -c(likobj$grad)*dlink(thetaTrans) +
        dpen(link(thetaTrans))*dlink(thetaTrans)
      likobj$info <- likobj$info*outer(dlink(thetaTrans),dlink(thetaTrans)) +
        ddpen(link(thetaTrans))*outer(dlink(thetaTrans),dlink(thetaTrans))
      return(likobj)
    }
    optObj <- fisher_scoring(objfun1, thetaTrans, link, T, 1e-3,
                             min(crtIter, maxIter - crtIter + 1))
    # 
    # optObj <- optim(theta, objfun1, dobjfun1, method = "L-BFGS-B",
    #                lower = c(0.01, rep(0, length(theta) - 2), 0.01^2), 
    #                control = list(maxit = min(crtIter, maxIter - crtIter + 1),
    #                               pgtol = 1e-3))
    # theta <- optObj$par
    crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
  }
  return(list(obj = optObj$loglik, parms = optObj$covparms))
}

varIdxLst <- lapply(c(1 : d), function(x){x})
criteriaObjLst <- lapply(varIdxLst, criteria)
objVec <- sapply(criteriaObjLst, function(x){x$obj})
bstVarIdx <- varIdxLst[[which.min(objVec)]] 
bstOptObj <- criteriaObjLst[[which.min(objVec)]]
while(length(bstVarIdx) < d)
{
  cat("Best variable indices are: ", bstVarIdx, "\n")
  idxPool <- setdiff(c(1 : d), bstVarIdx)
  varIdxLst <- lapply(idxPool, function(x){c(bstVarIdx, x)})
  criteriaObjLst <- lapply(varIdxLst, criteria)
  objVec <- sapply(criteriaObjLst, function(x){x$obj})
  bstOptObjNew <- criteriaObjLst[[which.min(objVec)]]
  if(stop_cond(bstOptObj, bstOptObjNew))
    break
  bstVarIdx <- varIdxLst[[which.min(objVec)]]
  bstOptObj <- bstOptObjNew
}
sink(file = NULL)








