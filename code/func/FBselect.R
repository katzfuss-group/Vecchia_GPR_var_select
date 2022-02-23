library(GpGp)
source("QCCD.R")

#' Forward-backward selection
#' 
#' @param idx crt selected covariates
#' @param theta crt values of cov parms
#' @param sr0 init sr
#' @param locsIn in-sample covariate matrix
#' @param locsOut out-sample covariate matrix
#' @param yIn in-sample obs vec
#' @param yOut out-sample obs vec
#' @param m Vecchia conditioning set size
#' @param k # selected covariates
#' @param lambda penalty parm
#' @param pen_fun penalty func
#' @param dpen_fun dpenalty func
#' @param ddpen_fun ddpenalty func
#' @param OOS_score out-of-sample score function
#' @param stop_con_fb forward-backward stop cond
#' @param miniQCCD mini-batch size for QCCD
#' @param miniGrad mini-batch size for gradient computation
#' @param convQCCD convergence level of QCCD
#' @param convCCD convergence level of CCD
#' @param cAmij Armijo constant
#' @param maxIterQCCD max # iterations of QCCD
#' @param maxIterCCD max # iterations of CCD
#' @param covFn covariance function name
#' @param minPosi a small positive number for lower bounds
#' @param mini bool for using mini-batching?
#' @param taper bool for using penalty tapering?
#' @param silent bool for silent execution?
fwd_bwd <- function(idx, theta, sr0, locsIn, locsOut, yIn, yOut, 
                    m, k, lambda, pen_fun, 
                    dpen_fun, ddpen_fun, OOS_score, stop_con_fb, miniQCCD, 
                    miniGrad, convQCCD = 1e-4, convCCD = 1e-4, cAmij = 1e-4, 
                    maxIterQCCD = 500, maxIterCCD = 40, 
                    covFn = "matern25_scaledim_sqrelevance", 
                    minPosi = 1e-8, mini = T, taper = F, silent = F, ...)
{
  d <- ncol(locsIn)
  # storage vars
  idxSet <- list()
  thetaSet <- list()
  scrVec <- c()
  iter <- 1
  # select kHat var each time
  while(length(idx) < d){
    if(!silent){
      cat("----------------------------------\n")
      cat("FB iter =", iter, "\n")
    }
    # gen lb for theta
    lb <- c(minPosi, rep(0, d), minPosi)
    lb[which.max(theta[2 : (d + 1)]) + 1] <- minPosi
    # MM and NN
    timePnt1 <- Sys.time()
    MMNNObj <- MM_NN(theta, idx, locsIn, yIn, lb, m)
    timePnt2 <- Sys.time()
    # QCCD
    QCCDObj <- QCCD_wrap(MMNNObj$thetaRel, MMNNObj$yOdr, MMNNObj$locsOdrRel, 
                         MMNNObj$NNarray, MMNNObj$lbRel, lambda, pen_fun, 
                         dpen_fun, ddpen_fun, convQCCD = convQCCD, 
                         convCCD = convCCD, cAmij = cAmij, 
                         maxIterQCCD = maxIterQCCD, maxIterCCD = maxIterCCD,
                         covFn = covFn, miniQCCD = miniQCCD, mini = mini, 
                         taper = taper, silent = silent, ...)
    timePnt3 <- Sys.time()
    theta <- get_theta(QCCDObj$covparms, idx, d)
    idx <- idx[sapply(idx, function(i){theta[i + 1] > 0})]
    idxSet[[iter]] <- idx
    thetaSet[[iter]] <- theta
    scrVec[iter] <- OOS_score(theta = theta, locsIn = locsIn, 
                             locsOut = locsOut, yIn = yIn, 
                             yOut = yOut, m = m, covFn = covFn)
    if(!silent){
      cat("FB iter =", iter, "\n")
      cat("idx:", idxSet[[iter]], "\n")
      cat("theta:", thetaSet[[iter]][c(1, idxSet[[iter]] + 1, d + 2)], "\n")
      cat("score:", scrVec[iter], "\n")
      cat("MM_NN used", as.numeric(timePnt2 - timePnt1, units = "secs"),
          "seconds\n")
      cat("QCCD used", as.numeric(timePnt3 - timePnt2, units = "secs"),
          "seconds\n")
    }
    # check stop
    iOpt <- stop_con_fb(scrVec, idxSet, thetaSet)
    if(iOpt > 0)
      break
    # compute grad
    if(mini){
      gradObj <- comp_grad_mean0(theta = theta, locs = MMNNObj$locsOdr, 
                                 y = MMNNObj$yOdr, 
                                 NNarray = MMNNObj$NNarray, miniGrad = miniGrad, 
                                 covFn = covFn)
    }else{
      gradObj <- comp_grad_mean0(theta = theta, locs = MMNNObj$locsOdr, 
                                 y = MMNNObj$yOdr, 
                                 NNarray = MMNNObj$NNarray, 
                                 miniGrad = length(MMNNObj$yOdr), 
                                 covFn = covFn)
    }
    # select var
    kHat <- min(k, d - length(idx))
    odrGrad <- order(gradObj$grad[2 : (d + 1)], decreasing = T)
    idxSel <- setdiff(odrGrad, idx)[1 : kHat]
    idx <- c(idx, idxSel)
    # init new var
    theta[idxSel + 1] <- sr0
    # next iter
    iter <- iter + 1
  }
  if(!silent){
    cat("----------------------------------\n")
  }
  if(iOpt < 1)
    iOpt <- length(scrVec)
  list(scr = scrVec[iOpt], idx = idxSet[[iOpt]], theta = thetaSet[[iOpt]])
}