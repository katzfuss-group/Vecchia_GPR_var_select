################################################################################
# Compare the bias of gradient estimates from four mini-batching methods
#     1. Randomly select 1 point and its batchsize - 1 nearest neighbors
#     2. Randomly select batchsize points
#     3. Randomly select batchsize rows in the NNarray
#     4. Stratified sampling batchsize rows in the NNarray
################################################################################

# devtools::install_github("https://github.com/SamCao1991/GpGp.git")
library(GpGp)
library(FNN)
library(lhs)
set.seed(123)
rm(list = ls())

# simulate GP
n <- 1e4
d <- 2
m <- 100
rSq <- c(1, 0.5)^2
sigmasq <- 1.0 # variance
tausq <- 0.05^2 # nugget
locs <- randomLHS(n, d)
locs <- locs * outer(rep(sqrt(n), n), 
                     1 / sqrt(colSums(locs^2)))
covM <- matern25_scaledim_sqrelevance(c(sigmasq, rSq, tausq), locs)
cholM <- t(chol(covM))
y <- as.vector(cholM %*% rnorm(n))
y <- y - mean(y)
theta <- c(sigmasq, rSq, tausq)
rm(covM, cholM)

# maximin order and NN conditioning
m <- min(m, n - 1)
locsScale <- locs %*% diag(sqrt(rSq), d, d)
odr <- order_maxmin(locsScale)
locsOdr <- locs[odr, , drop = F]
locsScaleOdr <- locsScale[odr, , drop = F]
yOdr <- y[odr]
NNarray <- find_ordered_nn(locsScaleOdr, m = m)

# batch size
batchSzVec <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 5000)
batchSzVec <- batchSzVec[batchSzVec <= n]

# Method 1
mtd1 <- function(batchSz){
  idx1 <- sample(x = 1 : n, size = 1)
  idx <- as.numeric(get.knnx(data = locsScale, 
                             query = locsScale[idx1, , drop = F],
                             k = batchSz)$nn.index)
  m <- min(m, batchSz - 1)
  locsBatch <- locs[idx, , drop = F]
  locsScaleBatch <- locsScale[idx, , drop = F]
  yBatch <- y[idx]
  odr <- order_maxmin(locsScaleBatch)
  locsBatchOdr <- locsBatch[odr, , drop = F]
  locsScaleBatchOdr <- locsScaleBatch[odr, , drop = F]
  yBatchOdr <- yBatch[odr]
  NNarray <- find_ordered_nn(locsScaleBatchOdr, m = m)
  
  vecchia_meanzero_loglik_grad_info(theta, 
                                    "matern25_scaledim_sqrelevance",
                                    yBatchOdr, locsBatchOdr, 
                                    NNarray)$grad
}

# Method 2
mtd2 <- function(batchSz){
  idx <- sample(x = 1 : n, size = batchSz, replace = F)
  m <- min(m, batchSz - 1)
  locsBatch <- locs[idx, , drop = F]
  locsScaleBatch <- locsScale[idx, , drop = F]
  yBatch <- y[idx]
  odr <- order_maxmin(locsScaleBatch)
  locsBatchOdr <- locsBatch[odr, , drop = F]
  locsScaleBatchOdr <- locsScaleBatch[odr, , drop = F]
  yBatchOdr <- yBatch[odr]
  NNarray <- find_ordered_nn(locsScaleBatchOdr, m = m)
  
  vecchia_meanzero_loglik_grad_info(theta, 
                                    "matern25_scaledim_sqrelevance",
                                    yBatchOdr, locsBatchOdr, 
                                    NNarray)$grad
}

# Method 3
mtd3 <- function(batchSz){
  idx <- sample(x = 1 : n, size = batchSz, replace = F)
  NNarrayBatch <- NNarray[idx, , drop = F]
  vecchia_meanzero_loglik_grad_info(theta, 
                                    "matern25_scaledim_sqrelevance",
                                    yOdr, locsOdr, 
                                    NNarrayBatch)$grad
}

# Method 4
mtd4 <- function(batchSz){
  idx <- sample(x = 1 : n, size = batchSz, replace = F, 
                prob = (1 : n)^(- 1 / d))
  NNarrayBatch <- NNarray[idx, , drop = F]
  vecchia_meanzero_loglik_grad_info(theta, 
                                    "matern25_scaledim_sqrelevance",
                                    yOdr, locsOdr, 
                                    NNarrayBatch)$grad
}

# Compute gradient with full sample
gradFull <- vecchia_meanzero_loglik_grad_info(theta, 
                                              "matern25_scaledim_sqrelevance",
                                              yOdr, locsOdr, 
                                              NNarray)$grad
covM <- matern25_scaledim_sqrelevance(c(sigmasq, rSq, tausq), locsOdr)
covMInv <- solve(covM)
dcovM <- d_matern25_scaledim_sqrelevance(theta, locsOdr)
gradExact <- rep(0, length(theta))
for(i in 1 : length(theta)){
  gradExact[i] <- t(yOdr) %*% covMInv %*% dcovM[, , i] %*% covMInv %*% yOdr -
    sum(covMInv * dcovM[, , i])
}
rm(covM, covMInv, dcovM)

# Compute bias for each method
nIter <- 1000
result <- array(dim = c(nIter, length(batchSzVec), (d + 2) * 2, 4), 
                dimnames = list(paste("Iter", 1 : nIter), 
                                paste("Batchsize", batchSzVec),
                                rep(c("Variance", paste("SR", 1 : d), "Nugget"), 
                                    2),
                                paste("Method", 1 : 4)))
for(j in 1 : length(batchSzVec)){
  batchSz <- batchSzVec[j]
  for(i in 1 : nIter){
    result[i, j, 1 : (d + 2), 1] <- mtd1(batchSz) / batchSz - gradFull / n
    result[i, j, 1 : (d + 2), 2] <- mtd2(batchSz) / batchSz - gradFull / n
    result[i, j, 1 : (d + 2), 3] <- mtd3(batchSz) / batchSz - gradFull / n
    result[i, j, 1 : (d + 2), 4] <- mtd4(batchSz) / batchSz - gradFull / n
    
    result[i, j, (d + 3) : ((d + 2) * 2), 1] <- mtd1(batchSz) / batchSz - 
      gradExact / n
    result[i, j, (d + 3) : ((d + 2) * 2), 2] <- mtd2(batchSz) / batchSz - 
      gradExact / n
    result[i, j, (d + 3) : ((d + 2) * 2), 3] <- mtd3(batchSz) / batchSz - 
      gradExact / n
    result[i, j, (d + 3) : ((d + 2) * 2), 4] <- mtd4(batchSz) / batchSz - 
      gradExact / n
  }
  cat("Batchsize = ", batchSz, " is done\n")
}

# save the results
save.image("grad_bias_four_minibatch.RData")

# plot the results
library(ggplot2)
library(scales)
for(i in 1 : ((d + 2) * 2)){
  df <- matrix(NA, length(batchSzVec) * 4, 4)
  varName <- dimnames(result)[[3]][i]
  varName <- gsub(" ", "", varName, fixed = T)
  for(j in 2 : 4){ # methods
    for(k in 1 : length(batchSzVec)){ # batch size
      df[(j - 1) * length(batchSzVec) + k, 1] <- batchSzVec[k]
      df[(j - 1) * length(batchSzVec) + k, 2] <- abs(mean(result[, k, i, j]))
      df[(j - 1) * length(batchSzVec) + k, 3] <-
        mean(result[, k, i, j]^2)
      df[(j - 1) * length(batchSzVec) + k, 4] <- j
    }
  }
  colnames(df) <- c("batchsz", "bias", "MSE", "mtdID")
  df <- as.data.frame(df)
  df$mtdID <- paste("Method", df$mtdID)
  df <- df[complete.cases(df), ]
  plt <- ggplot(data = df, aes(x = batchsz, y = bias, col = mtdID, 
                               lty = mtdID)) +
    geom_line() +
    scale_x_continuous(trans = pseudo_log_trans(sigma = 1),
                       breaks = batchSzVec) +
    scale_color_manual(breaks = unique(df$mtdID),
                       values = rainbow(4)) +
    theme(legend.position = "none")
  if(i != 1 && i != d + 2 && i != d + 3 && i != 2*d + 4)
    plt <- plt + scale_y_continuous(trans = pseudo_log_trans(sigma = 0.002))
  if(i <= d + 2){
    ggsave(paste0("grad_bias_", varName, ".pdf"), plot = plt,
           width = 7, height = 5)
  }else{
    ggsave(paste0("grad_bias_exact_", varName, ".pdf"), plot = plt,
           width = 7, height = 5)
  }
  

  plt <- ggplot(data = df, aes(x = batchsz, y = MSE, col = mtdID, 
                               lty = mtdID)) +
    geom_line() +
    scale_x_continuous(trans = pseudo_log_trans(sigma = 1),
                       breaks = batchSzVec) +
    scale_color_manual(breaks = unique(df$mtdID),
                       values = rainbow(4)) + 
    theme(legend.position = "none")
  if(i != 1 && i != d + 2 && i != d + 3 && i != 2*d + 4)
    plt <- plt + scale_y_continuous(trans = pseudo_log_trans(sigma = 0.0001))
  if(i <= d + 2){
    ggsave(paste0("grad_MSE_", varName, ".pdf"), plot = plt,
           width = 7, height = 5)
  }else{
    ggsave(paste0("grad_MSE_exact_", varName, ".pdf"), plot = plt,
           width = 7, height = 5)
  }
}



























