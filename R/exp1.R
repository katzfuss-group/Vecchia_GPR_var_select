# Experiment 1

n <- 1e4
d <- 20
r <- c(10, 5, 2, 1, 0.5, rep(0, d - 5))
sigma <- 1.0 # variance
tausq <- 0.05^2 # nugget

set.seed(123)
locs <- lhs::randomLHS(n, d)
covM <- GpGp::matern25_scaledim(c(sigma, 1 / r, tausq), locs)