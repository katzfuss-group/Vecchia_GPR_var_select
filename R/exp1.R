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
# For 'matern_scaledim'
# nuInit <- 1.5
# theta <- c(sigmasqInit, rInit, nuInit, tausqInit)
theta <- c(sigmasqInit, rInit, tausqInit)
crtIter <- 1
maxIter <- 100
m <- 30
lambda <- 1
while(maxIter >= crtIter)
{
  locsScal <- locs %*% diag(theta[2 : (d + 1)])
  odr <- GpGp::order_maxmin(locsScal)
  yOdr <- y[odr]
  locsOdr <- locs[odr, ]
  XOdr <- X[odr, , drop = F]
  NNarray <- GpGp::find_ordered_nn(locsOdr, m = m)

  objfun1 <- function(theta){
    GpGp::vecchia_profbeta_loglik_grad_info(theta, "matern25_scaledim_relevance",
                                            yOdr, XOdr, locsOdr, 
                                            NNarray)
  }
  theta <- quad_cdsc_L1(objfun1, theta, lambda, c(2 : (d + 1)), 1e-3,
               max_iter = min(crtIter, maxIter - crtIter + 1))$covparms
  crtIter <- crtIter * 2
}
















