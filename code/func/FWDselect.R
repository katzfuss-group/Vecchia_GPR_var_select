library(parallel)
source("fisher_scoring.R")
source("misc.R")

#' Regularization path using Vecc log-lik and some penalty
#' 
#' @param var0 init var
#' @param tau0 init tau
#' @param sr0 init sr
#' @param locs covariate matrix
#' @param y obs vec
#' @param m Vecchia conditioning set size
#' @param OOS_score out-of-sample score function
#' @param stop_con_path path stop cond
#' @param pIn in-sample proportion
#' @param spltSeed locs and y random split seed
#' @param conv convergence level of Fisher scoring
#' @param maxIter max # iterations of Fisher scoring
#' @param covFn covariance function name
#' @param link_func link func, e.g., log
#' @param resp_func response func, e.g., exp
#' @param dresp_func dresponse func, e.g., exp
#' @param silent bool for silent execution?
forward_selection_fisher <- function(var0, tau0, sr0, locs, y, m, OOS_score, 
                                     stop_con_path, pIn = 0.25, spltSeed = NULL, 
                                     conv = 1e-4, maxIter = 100, 
                                     covFn = "matern25_scaledim_sqrelevance",
                                     link_func = log, resp_func = exp, 
                                     dresp_func = exp,
                                     silent = F, cluster = NULL)
{
  n <- nrow(locs)
  d <- ncol(locs)
  
  # split locs and y
  spltObj <- splt_locs_y(locs = locs, y = y, p1 = pIn, seed = spltSeed)
  locsIn <- spltObj$locs1
  yIn <- spltObj$y1
  locsOut <- spltObj$locs2
  yOut <- spltObj$y2
  
  idxSet <- list()
  thetaSet <- list()
  scrVec <- c()
  selectIdx <- c()
  
  theta <- c(var0, rep(sr0, d), tau0)
  for(i in 1 : d)
  {
    if(!silent){
      cat("======================================================\n")
      cat("i =", i, "\n")
    }
    bgnTime <- Sys.time()
    unselectIdx <- setdiff(1 : d, selectIdx)
    idxList <- lapply(unselectIdx, function(x){c(selectIdx, x)})
    fisher_scoring_meanzero_wrap_par <- function(idx, ...){
      source("../func/misc.R", chdir = T)
      source("../func/fisher_scoring.R", chdir = T)
      fisher_scoring_meanzero_wrap(idx, ...)
    }
    clusterExport(cluster, ls(), envir = environment())
    models <- parLapply(cluster, idxList, fisher_scoring_meanzero_wrap_par, 
                        m = m, theta = theta, locs = locsIn, y = yIn, 
                        covFn = covFn, maxIter = maxIter, conv = conv, 
                        link_func = link_func, resp_func = resp_func, 
                        dresp_func = dresp_func, silent = silent)
    scores <- parLapply(cluster, models, OOS_score, 
                        locsIn = locsIn, locsOut = locsOut,
                        yIn = yIn, yOut = yOut, m = m, covFn = covFn)
    scores <- unlist(scores)
    idxBest <- which.min(scores)
    endTime <- Sys.time()
    
    # Store results
    idxSet[[i]] <- idxList[[idxBest]]
    thetaSet[[i]] <- models[[idxBest]]
    scrVec[i] <- scores[idxBest]
    theta[c(1, idxSet[[i]] + 1, d + 2)] <- 
      models[[idxBest]][c(1, idxSet[[i]] + 1, d + 2)]
    selectIdx <- idxList[[idxBest]]
    
    if(!silent){
      cat("i =", i, ":\n")
      cat("idx:", idxSet[[i]], "\n")
      cat("theta:", thetaSet[[i]][c(1, idxSet[[i]] + 1, d + 2)], "\n")
      cat("score:", scrVec[i], "\n")
      cat("Current iteration used", as.numeric(endTime - bgnTime, 
                                               units = "secs"), 
          "seconds\n")
    }
    
    idxFWD <- stop_con_path(scrVec, idxSet, thetaSet)
    if(idxFWD > 0)
      break
  }
  if(!silent){
    cat("======================================================\n")
    cat("Results of forward selection:\n")
    cat("idx:", idxSet[[idxFWD]], "\n")
    cat("theta:", thetaSet[[idxFWD]][c(1, idxSet[[idxFWD]] + 1, d + 2)], "\n")
    cat("score:", scrVec[idxFWD], "\n")
  }
  list(scr = scrVec[idxFWD], idx = idxSet[[idxFWD]], 
       theta = thetaSet[[idxFWD]])
}