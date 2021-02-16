#' Batch-start the optimization over the relevance parameters 
#' 
#' @param nonrel_parms init val for non-relevance parms
#' @param d number of relevance parms
#' @param batch_sz batch size
#' @param extra_batch extra iterations for randomly selecting batch_sz zero relevance parms and add one
#' @param opt_fun the optimization function
#' @param ... Other arguments to opt_fun
batch_start <- function(nonrel_parms, d, batch_sz, extra_batch, opt_fun, rel_parms = rep(0, d), ...)
{
  step <- floor(d / batch_sz)
  parms <- c(nonrel_parms[1], rel_parms, nonrel_parms[-1])
  if(step > 0)
  {
    for(i in 1 : step)
    {
      idx <- seq(from = i, to = d, by = step)
      parms[idx + 1] <- parms[idx + 1] + 1
      optObj <- opt_fun(parms, ...)
      parms <- optObj$covparms
      cat("step", i, "non-zero parm indices are", which(parms[2 : (d + 1)] > 0), "\n")
    }
    # for(i in 1 : extra_batch)
    # {
    #   
    # }
  }else{
    parms <- c(nonrel_parms[1], rel_parms + 1, nonrel_parms[-1])
    optObj <- opt_fun(parms, ...)
    parms <- optObj$covparms
    cat("step", 1, "non-zero parm indices are", which(parms[2 : (d + 1)] > 0), "\n")
  }
  
  optObj <- opt_fun(parms, ...)
  cat("After a final optimization, non-zero parm indices are", which(parms[2 : (d + 1)] > 0), "\n")
  return(optObj)
}