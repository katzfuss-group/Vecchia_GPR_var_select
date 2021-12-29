

#' REG path stop condition based on 1% impr OOS score
#' 
#' @param scrVec score vector
#' @param idxSet idx list
#' @param thetaSet theta list
#' @param minPosi a small positive number checking if theta is meaningful
stop_con_OOS1_path <- function(scrVec, idxSet, thetaSet, minPosi = 1e-4)
{
  i <- length(scrVec)
  if(i < 2)
    return(-1)
  d <- length(thetaSet[[1]]) - 2
  if(scrVec[i] > scrVec[i - 1] * 0.99 &
     sum(thetaSet[[i - 1]][2 : (d + 1)]) > minPosi &
     sum(thetaSet[[i]][2 : (d + 1)]) > minPosi
     )
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
  if(scrVec[i] > scrVec[i - 1] * 0.99)
    return(i - 1)
  else
    return(-1)
}