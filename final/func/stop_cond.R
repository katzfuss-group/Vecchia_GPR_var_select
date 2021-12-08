

#' REG path stop condition based on 1% impr OOS score
#' 
#' @param scrVec score vector
#' @param idxSet idx list
#' @param thetaSet theta list
stop_con_OOS1_path <- function(scrVec, idxSet, thetaSet)
{
  i <- length(scrVec)
  if(i < 2)
    return(-1)
  d <- length(thetaSet[[1]]) - 2
  if(scrVec[i] > scrVec[i - 1] * 0.99 &
     sum(thetaSet[[i - 1]][2 : (d + 1)]) > 0.1 &
     sum(thetaSet[[i]][2 : (d + 1)]) > 0.1)
    return(i - 1)
  else
    return(-1)
}

#' Forward-backward stop condition based on 1% impr OOS score
#' 
#' @param scrVec score vector
#' @param idxSet idx list
#' @param thetaSet theta list
stop_con_OOS1_fb <- function(scrVec, idxSet, thetaSet)
{
  i <- length(scrVec)
  if(i < 2)
    return(-1)
  # from the 2nd iter, if SRs are still too small, stop crt fwd bwd selection
  # meaning that lambda is too big
  if(sum(thetaSet[[i - 1]][2 : (d + 1)]) < 0.1)
    return(i)
  d <- length(thetaSet[[1]]) - 2
  if(scrVec[i] > scrVec[i - 1] * 0.99 &
     sum(thetaSet[[i - 1]][2 : (d + 1)]) > 0.1 &
     sum(thetaSet[[i]][2 : (d + 1)]) > 0.1)
    return(i - 1)
  else
    return(-1)
}