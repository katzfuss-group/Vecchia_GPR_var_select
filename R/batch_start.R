for(i in 1 : 5)
{
  idx <- i + seq(1, d, by = d / 20)
  theta[idx] <- theta[idx] + 1
  theta <- quad_cdsc_L1(objfun, objfun_gdfm, locsOdr, 1, theta, lambda, 1e-3, silent = T, 
                        max_iter = innerIter, max_iter2 = 40, lb_nonrel_parms = lb_nonrel_parms)$covparms
  cat(i, ": estimated theta = ", theta, "\n")
}


thetaTrans <- fisher_scoring(objfun1, thetaTrans, link, T, 1e-3,
                             100)$logparms

batch_start <- function(nonrel_parms, d, batch_sz, opt_fun, ...)
{
  
}