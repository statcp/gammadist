#' Construct lower confidence limit, upper confidence limit and confidence interval for the
#' parameters of the gamma distribution.
#'
#' @param x observations of a gamma distribution.
#' @param alpha value of alpha, such that an (1-alpha) lower confidence limit, an (1-alpha) upper confidence
#' limit or an (1-alpha) confidence interval is returned. If alpha is not specified, the default
#' value is 0.05.
#' @param B number of realizations of the GPQs. If B is not specified, the default value is 2000.
#'
#' @details Assume that x_1, ..., x_n from a gamma distribution are available. 
#' Consider two statistics
#' L = L(x_1, ..., x_n) and U = U(x_1, ..., x_n). If
#' P(L <= k <= U) = 1-alpha,
#' then we call the interval [L, U] an (1-alpha) confidence interval for the shape parameter k. When U
#' equals infinity, we call L an (1-alpha) lower confidence limit for the shape parameter. Similarly,
#' when L equals -infnity, we call U an (1-alpha) upper confidence limit for the shape parameter.
#' The same arguments hold for the scale and rate parameter.
#'
#' @return conflimits returns a dataframe which contains a list of the variables shape, scale and rate. 
#' Each variable is a vector containing 4 values, where the first two are the equal-tailed 1-alpha confidence intervals, 
#' and the last two are the 1-alpha lower and upper confidence limits, respectively.
#' @export
#'
#' @references Wang, B. X. and Wu, F. (2018), “Inference on the Gamma Distribution,” Technometrics, 60(2), 235–244.
#' 
#' @examples x <- rGamma(100, shape = 3, rate = 2)
#' conf.all <- conflimits(x, alpha=0.1, B=5000)
#' shape.conf <- conf.all$shape



conflimits <- function(x,alpha=0.05,B=2000) {
  gpq.all <- pargpq (x,B)
  conf.all <- apply(gpq.all,2,function(x) quantile(x,c(alpha/2,1-alpha/2,alpha,1-alpha)))
  conf.all <- as.data.frame(conf.all,row.names=c('low-int','up-int','low-lim','up-lim'))
  return(conf.all)
}
