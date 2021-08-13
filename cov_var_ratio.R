set.seed(123)

n <- 5e6
rsqVec <- c((1 : 100) / 10)^2
ratioVec <- rep(NA, length(rsqVec))
X <- matrix(rnorm(n * 3), n, 3)
ratioThres <- cov((X[, 1] - X[, 2])^2, (X[, 1] - X[, 3])^2) / 
  var((X[, 1] - X[, 2])^2)
for(i in 1 : length(rsqVec))
{
  rsq <- rsqVec[i]
  ratioVec[i] <- cov(- (X[, 1] - X[, 2])^2, exp(- rsq * (X[, 1] - X[, 3])^2)) / 
    cov(- (X[, 1] - X[, 2])^2, exp(- rsq * (X[, 1] - X[, 2])^2))
}
ratioVec

library(ggplot2)
mydf <- data.frame(rsq = rsqVec, ratio = ratioVec)
ggplot(data = mydf, aes(x = rsq, y = ratio)) +
  geom_line() +
  geom_hline(yintercept = ratioThres, size = 0.5, 
             color = "red", lty = "dashed") 
ggsave(paste0("normal_locs_ratio.pdf"), width = 7, height = 5)




set.seed(123)

n <- 5e6
rsqVec <- c((1 : 100) / 10)^2
ratioVec <- rep(NA, length(rsqVec))
X <- matrix(runif(n * 3, min = -0.5, max = 0.5), n, 3)
ratioThres <- cov((X[, 1] - X[, 2])^2, (X[, 1] - X[, 3])^2) / 
  var((X[, 1] - X[, 2])^2)
for(i in 1 : length(rsqVec))
{
  rsq <- rsqVec[i]
  ratioVec[i] <- cov(- (X[, 1] - X[, 2])^2, exp(- rsq * (X[, 1] - X[, 3])^2)) / 
    cov(- (X[, 1] - X[, 2])^2, exp(- rsq * (X[, 1] - X[, 2])^2))
}
ratioVec

library(ggplot2)
mydf <- data.frame(rsq = rsqVec, ratio = ratioVec)
ggplot(data = mydf, aes(x = rsq, y = ratio)) +
  geom_line() +
  geom_hline(yintercept = ratioThres, size = 0.5, 
             color = "red", lty = "dashed") 
ggsave(paste0("uniform_locs_ratio.pdf"), width = 7, height = 5)


