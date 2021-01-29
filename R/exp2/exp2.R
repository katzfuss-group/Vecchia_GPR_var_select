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
lambdaVec <- exp(seq(0, log(n), length.out = 5))
lb_nonrel_parms <- c(0.01^2, 0.01^2)

# quad_cdsc_L1 method
theta <- c(sigmasqInit, rInit, tausqInit)
innerIter <- 20
outerIter <- 5
m <- 30
sink("quad_cdsc_L1.out")
for(i in 1 : outerIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  # Define functions for cross-validation
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
                 max_iter = innerIter, max_iter2 = 40, lb_nonrel_parms = lb_nonrel_parms)
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
  
  loss <- rep(NA, length(lambdaVec))
  idxRnd <- sample(c(1 : n), n, F)
  for(j in 1 : length(lambdaVec))
  {
    loss[j] <- cross_valid(est_func, pred_func, crit_MSE, locsOdr, XOdr, yOdr, 
                           m, lambdaVec[j], 5, idxRnd)
    cat(lambdaVec[j], ": loss ", loss[j], "\n")
  }
  lambda <- lambdaVec[which.min(loss)]
  cat(i, ": best lambda = ", lambda, "\n")
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
  
  theta <- quad_cdsc_L1(objfun, objfun_gdfm, locsOdr, 1, theta, lambda, 1e-3, silent = F, 
               max_iter = innerIter, max_iter2 = 40, lb_nonrel_parms = lb_nonrel_parms)$covparms
  cat(i, ": estimated theta = ", theta, "\n")
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
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(link(thetaTrans[2 : (d + 1)]))
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)

  objfun1 <- function(thetaTrans){
    cat("c")
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
  thetaTrans <- fisher_scoring(objfun1, thetaTrans, link, T, 1e-3,
                        min(crtIter, maxIter - crtIter + 1))$logparms
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# optim package from R
theta <- c(sigmasqInit, rInit, tausqInit)
crtIter <- 1
maxIter <- 100
m <- 30
sink("L-BFGS-B.out")
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  
  objfun1 <- function(theta){
    cat(theta, "\n")
    likobj <- GpGp::vecchia_profbeta_loglik(theta, 
                                            "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
    return(-likobj$loglik + lambda * sum(theta[2 : (d + 1)]))
  }
  dobjfun1 <- function(theta){
    # cat("c")
    cat(theta, "\n")
    likobj <- GpGp::vecchia_profbeta_loglik_grad_info(theta, 
                                                      "matern25_scaledim_relevance",
                                                      yOdr, XOdr, locsOdr, 
                                                      NNarray)
    return(-likobj$grad + lambda * c(0, rep(1, d), 
                                     rep(0, length(theta) - d - 1)))
  }
  
  theta <- optim(theta, objfun1, dobjfun1, method = "L-BFGS-B",
                 lower = c(0.01, rep(0, d), 0.01^2), 
                 control = list(maxit = min(crtIter, maxIter - crtIter + 1),
                                trace = 0,
                                pgtol = 1e-3))$par
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# Forward selection
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
      cat("c")
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
  return(list(obj = optObj$value, parms = optObj$par))
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








