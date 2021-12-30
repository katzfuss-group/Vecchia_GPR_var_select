library(parallel)
library(Rcgmin)
library(GpGp)
library(mvtnorm)
source("penalty.R")
source("misc.R")

#' Regularization path using Yi(2011) method
#' @param theta_gen theta init generator
#' @param ngen number of random starts
#' @param locs covariate matrix
#' @param y obs vec
#' @param m Vecchia conditioning set size
#' @param lambVec penalty parms vec
#' @param pen_fun penalty func
#' @param dpen_fun dpenalty func
#' @param ddpen_fun ddpenalty func
#' @param OOS_score out-of-sample score function
#' @param stop_con_path path stop cond
#' @param pIn in-sample proportion
#' @param spltSeed locs and y random split seed
#' @param covFn covariance function name
#' @param maxIter max # iterations of Fisher scoring
#' @param conv convergence level of Fisher scoring
#' @param link_func link func, e.g., log
#' @param resp_func response func, e.g., exp
#' @param dresp_func dresponse func, e.g., exp
#' @param silent bool for silent execution?
PGPR_Yi11 <- function(theta_gen, ngen, locs, y, m, lambVec, pen_fun, 
                      dpen_fun, ddpen_fun, OOS_score, stop_con_path, pIn, 
                      spltSeed = NULL, covFn = "matern25_scaledim_sqrelevance", 
                      maxIter = 100, conv = 1e-4, 
                      link_func = log, resp_func = exp, dresp_func = exp, 
                      silent = T)
{
    d <- ncol(locs)
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
    # obj funcs
    obj_func <- function(theta)
    {
        covM <- get(covFn)(resp_func(theta), locsIn)
        - dmvnorm(yIn, sigma = covM, log = T) + pen_fun(resp_func(theta), lambda)
    }
    grad_func <- function(theta)
    {
        covM <- get(covFn)(resp_func(theta), locsIn)
        dcovM <- get(paste0("d_", covFn))(resp_func(theta), locsIn)
        covMInv <- solve(covM)
        covMInvy <- covMInv %*% yIn
        idx <- 1 : length(theta)
        grad <- 
            unlist(mclapply(idx, 
                            function(x){
                                (sum(covMInv * dcovM[, , x]) - 
                                     sum(t(covMInvy) %*% dcovM[, , x] %*% 
                                             covMInvy)) / 2},
                            mc.cores = detectCores() - 1))
        grad * dresp_func(theta) + dpen_fun(resp_func(theta), lambda) * 
            dresp_func(theta)
    }
    for(i in 1 : length(lambVec))
    {
        lambda = lambVec[i]
        if(!silent){
            cat("======================================================\n")
            cat("lambda =", lambda, "\n")
        }
        # loop over different init val
        for(j in 1 : ngen){
            if(!silent){
                cat("CG with the", j, "th initial value\n")
            }
            theta <- theta_gen()
            optObj <- Rcgmin(par = link_func(theta), fn = obj_func, 
                             gr = grad_func, 
                             control = list(maxit = maxIter, trace = !silent))
            if(j == 1){
                thetaFit <- optObj$par
                obj <- optObj$value
            }else{
                if(optObj$value < obj){
                    thetaFit <- optObj$par
                    obj <- optObj$value
                }
            }
        }
        # Store results
        idxSet[[i]] <- which(resp_func(thetaFit[2 : (d + 1)]) > 1e-7)
        thetaSet[[i]] <- resp_func(thetaFit) 
        scrVec[i] <-  OOS_score(theta = resp_func(thetaFit), locsIn = locsIn, 
                                locsOut = locsOut, yIn = yIn, 
                                yOut = yOut, m = m, covFn = covFn)
        # check stop
        iOpt <- stop_con_path(scrVec, idxSet, thetaSet)
        if(iOpt > 0)
            break
    }
    if(iOpt < 1)
        iOpt <- length(scrVec)
    if(!silent){
        cat("======================================================\n")
        cat("Results of PGPR (Yi2011):\n")
        cat("lambda:", lambVec[iOpt], "\n")
        cat("idx:", idxSet[[iOpt]], "\n")
        cat("theta:", thetaSet[[iOpt]][c(1, idxSet[[iOpt]] + 1, d + 2)], "\n")
        cat("score:", scrVec[iOpt], "\n")
    }
    list(lambda = lambVec[iOpt], scr = scrVec[iOpt], idx = idxSet[[iOpt]], 
         theta = thetaSet[[iOpt]])
}











