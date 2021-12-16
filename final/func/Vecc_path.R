source("FBselect.R")

#' Regularization path using Vecc log-lik and some penalty
#' 
#' @param var0 init var
#' @param tau0 init tau
#' @param sr0 init sr
#' @param locs covariate matrix
#' @param y obs vec
#' @param m Vecchia conditioning set size
#' @param k # selected covariates
#' @param lambVec penalty parms vec
#' @param pen_fun penalty func
#' @param dpen_fun dpenalty func
#' @param ddpen_fun ddpenalty func
#' @param OOS_score out-of-sample score function
#' @param stop_con_path path stop cond
#' @param stop_con_fb forward-backward stop cond
#' @param miniCQCD mini-batch size for CQCD
#' @param miniGrad mini-batch size for gradient computation
#' @param pIn in-sample proportion
#' @param spltSeed locs and y random split seed
#' @param convCQCD convergence level of CQCD
#' @param convCCD convergence level of CCD
#' @param cAmij Armijo constant
#' @param maxIterCQCD max # iterations of CQCD
#' @param maxIterCCD max # iterations of CCD
#' @param covFn covariance function name
#' @param minPosi a small positive number for lower bounds
#' @param mini bool for using mini-batching?
#' @param taper bool for using penalty tapering?
#' @param silent bool for silent execution?
Vecc_path <- function(var0, tau0, sr0, locs, y, m, k, lambVec, pen_fun, 
                      dpen_fun, ddpen_fun, OOS_score, stop_con_path,
                      stop_con_fb, miniCQCD, miniGrad, pIn = 0.25, 
                      spltSeed = NULL, convCQCD = 1e-4, 
                      convCCD = 1e-4, cAmij = 1e-4, maxIterCQCD = 500, 
                      maxIterCCD = 40, covFn = "matern25_scaledim_sqrelevance", 
                      minPosi = 1e-8, mini = T, taper = F, silent = F, ...)
{
  n <- nrow(locs)
  d <- ncol(locs)
  # init theta and idx
  rndIdx <- sample(1 : d, 1)
  theta <- c(var0, rep(0, d), tau0)
  theta[rndIdx + 1] <- sr0
  idx <- rndIdx
  # storage vars
  idxSet <- list()
  thetaSet <- list()
  scrVec <- c()
  # split locs and y
  spltObj <- splt_locs_y(locs = locs, y = y, p1 = pIn, seed = spltSeed)
  locsIn <- spltObj$locs1
  yIn <- spltObj$y1
  locsOut <- spltObj$locs2
  yOut <- spltObj$y2
  # loop over lambda
  for(i in 1 : length(lambVec)){
    lambda <- lambVec[i]
    if(!silent){
      cat("======================================================\n")
      cat("lambda =", lambda, "\n")
    }
    # re-init var and tau for each lambda
    theta[1] <- var0
    theta[d + 2] <- tau0
    # 1st level opt
    bgnTime <- Sys.time()
    optObj <- fwd_bwd(idx = idx, theta = theta, sr0 = sr0, locsIn = locsIn, 
                      locsOut = locsOut, yIn = yIn, yOut = yOut, m = m, k = k, 
                      lambda = lambda, pen_fun = pen_fun, dpen_fun = dpen_fun, 
                      ddpen_fun = ddpen_fun, OOS_score = OOS_score, 
                      stop_con_fb = stop_con_fb, miniCQCD = miniCQCD, 
                      miniGrad = miniGrad, convCQCD = convCQCD, 
                      convCCD = convCCD, cAmij = cAmij, 
                      maxIterCQCD = maxIterCQCD, maxIterCCD = maxIterCCD, 
                      covFn = covFn, minPosi = minPosi, mini = mini, 
                      taper = taper, silent = silent, ...)
    endTime <- Sys.time()
    # save results
    idxSet[[i]] <- optObj$idx
    thetaSet[[i]] <- optObj$theta
    scrVec[i] <- optObj$scr
    idx <- optObj$idx
    theta <- optObj$theta
    if(!silent){
      cat("Results for lambda =", lambVec[i], ":\n")
      cat("idx:", idxSet[[i]], "\n")
      cat("theta:", thetaSet[[i]][c(1, idxSet[[i]] + 1, d + 2)], "\n")
      cat("score:", scrVec[i], "\n")
      cat("fwd_bwd used", as.numeric(endTime - bgnTime, units = "secs"), 
          "seconds\n")
    }
    # check stop
    iOpt <- stop_con_path(scrVec, idxSet, thetaSet)
    if(iOpt > 0)
      break
  }
  if(iOpt < 1)
    iOpt <- length(scrVec)
  if(!silent){
    cat("======================================================\n")
    cat("Results of Vecchia regularization path:\n")
    cat("lambda:", lambVec[iOpt], "\n")
    cat("idx:", idxSet[[iOpt]], "\n")
    cat("theta:", thetaSet[[iOpt]][c(1, idxSet[[iOpt]] + 1, d + 2)], "\n")
    cat("score:", scrVec[iOpt], "\n")
  }
  list(lambda = lambVec[iOpt], scr = scrVec[iOpt], idx = idxSet[[iOpt]], 
       theta = thetaSet[[iOpt]])
}