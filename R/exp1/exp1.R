# Experiment 1
n <- 1e4
d <- 20
r <- c(10, 5, 2, 1, 0.5, rep(0, d - 5))
sigmasq <- 1.0 # variance
tausq <- 0.05^2 # nugget

set.seed(123)
locs <- lhs::randomLHS(n, d)
covM <- GpGp::matern25_scaledim_relevance(c(sigmasq, r, tausq), locs)
cholM <- t(chol(covM))
y <- as.vector(cholM %*% rnorm(n))
X <- matrix(1, n, 1)

rInit <- rep(1, d)
sigmasqInit <- 0.25
tausqInit <- 0
lambda <- 1

# quad_cdsc_L1 method
theta <- c(sigmasqInit, rInit, tausqInit)
crtIter <- 1
maxIter <- 100
m <- 30
sink("quad_cdsc_L1.out")
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)

  objfun1 <- function(theta){
    GpGp::vecchia_profbeta_loglik(theta, "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
  }
  objfun1GDFM <- function(theta){
    cat("c")
    GpGp::vecchia_profbeta_loglik_grad_info(theta, "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
  }
  
  theta <- quad_cdsc_L1(objfun1, objfun1GDFM, theta, lambda, c(2 : (d + 1)), 1e-3, silent = T, 
               max_iter = min(crtIter, maxIter - crtIter + 1), max_iter2 = 40)$covparms
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}
sink(file = NULL)

# Fisher scoring
tausqInitFS <- tausqInit + 0.01^2
theta <- c(sigmasqInit, rInit, tausqInitFS)
crtIter <- 1
maxIter <- 100
m <- 30

linkfuns <- get_linkfun("matern25_scaledim")
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

# Fisher scoring with range parameters
tausqInitFS <- tausqInit + 0.01^2
theta <- c(sigmasqInit, 1 / rInit, tausqInitFS)
crtIter <- 1
maxIter <- 100
m <- 30

linkfuns <- get_linkfun("matern25_scaledim")
link <- linkfuns$link
dlink <- linkfuns$dlink
invlink <- linkfuns$invlink
thetaTrans <- invlink(theta)

pen <- function(theta){lambda * sum(1 / theta[2 : (d + 1)])}
dpen <- function(theta){lambda * c(0, 
                                   -lambda / theta[2 : (d + 1)]^2, 
                                   rep(0, length(theta) - d - 1))}
ddpen <- function(theta){diag(c(0, 
                                2 * lambda / theta[2 : (d + 1)]^3,
                                rep(0, length(theta) - d - 1)))}
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(1 / link(thetaTrans[2 : (d + 1)]))
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)
  
  objfun1 <- function(thetaTrans){
    likobj <- 
      GpGp::vecchia_profbeta_loglik_grad_info(link(thetaTrans), 
                                              "matern25_scaledim",
                                              yOdr, XOdr, locsOdr, 
                                              NNarray)
    likobj$loglik <- -likobj$loglik + pen(link(thetaTrans))
    likobj$grad <- -c(likobj$grad)*dlink(thetaTrans) +
      dpen(link(thetaTrans))*dlink(thetaTrans)
    likobj$info <- likobj$info*outer(dlink(thetaTrans),dlink(thetaTrans)) +
      ddpen(link(thetaTrans))*outer(dlink(thetaTrans),dlink(thetaTrans))
    return(likobj)
  }
  thetaTrans <- fisher_scoring(objfun1, thetaTrans, link, F, 1e-3,
                               min(crtIter, maxIter - crtIter + 1))$logparms
  crtIter <- crtIter + min(crtIter, maxIter - crtIter + 1)
}

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
    cat("c")
    likobj <- GpGp::vecchia_profbeta_loglik_grad_info(theta, 
                                                      "matern25_scaledim_relevance",
                                                      yOdr, XOdr, locsOdr, 
                                                      NNarray)
    return(-likobj$loglik + lambda * sum(theta[2 : (d + 1)]))
  }
  dobjfun1 <- function(theta){
    likobj <- GpGp::vecchia_profbeta_loglik_grad_info(theta, 
                                                      "matern25_scaledim_relevance",
                                                      yOdr, XOdr, locsOdr, 
                                                      NNarray)
    return(-likobj$grad + lambda * c(0, rep(1, d), 
                                     rep(0, length(theta) - d - 1)))
  }
  
  theta <- optim(theta, objfun1, dobjfun1, method = "L-BFGS-B",
                 lower = c(0.01, rep(0, length(theta) - 1)), 
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








