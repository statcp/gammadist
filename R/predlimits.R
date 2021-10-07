#' Construct lower prediction limit, upper prediction limit and prediction interval for a
#' single future measurement of the gamma distribution.
#'
#' @param x observations of a gamma distribution.
#' @param alpha value of alpha, such that an (1-alpha) lower prediction limit, an (1-alpha) upper prediction
#' limit or an (1-alpha) prediction interval is returned. If alpha is not specified, the default
#' value is 0.05.
#' @param B number of realizations of the GPQs. If B is not specified, the default value is 2000.
#'
#' @details Assume that X_1, ..., X_n from a gamma distribution are available. 
#' Consider two statistics PL = PL(X_1, ..., X_n) and PU = PU(X_1, ..., X_n) and a future measurement
#' X_(n+1) from the same gamma distribution. If PL, PU and X_(n+1) have the following relationship,
#' P(PL <= X_(n+1) <= PU) = 1-alpha,
#' where 0 < alpha < 1 is a constant, then [PL,PU] is an (1-alpha) prediction interval for a single future measurement.
#' When PU equals infinity, we call PL an (1-alpha) lower prediction limit for a single future measurement.
#' Similarly, when PL equals -infinity, we call PU an (1-alpha) upper prediction limit for a single future measurement.
#'
#' @return predlimits returns a dataframe containing one vector where the first two elements make for the (1-alpha) prediction interval, the third element is the 
#' (1-alpha) lower prediction limit and the last is the (1-alpha) upper prediction limit.
#' @export
#'
#' @references Wang, B. X. and Wu, F. (2018), “Inference on the Gamma Distribution,” Technometrics, 60(2), 235–244.
#' @examples x <- rGamma(100, shape = 4, rate = 3)
#' predlimits(x, alpha=0.1, B=5000)

predlimits <- function(x,alpha=0.05,B=2000) {
  gpq.all <- pargpq (x,B)
  pred.all <- mapply(function(shape,scale) rGamma(1,shape=shape,scale=scale), shape=gpq.all$shape,scale=gpq.all$scale)
  pred <- quantile(pred.all,c(alpha/2,1-alpha/2,alpha,1-alpha))
  return(data.frame(pred,row.names=c('low-int','up-int','low-lim','up-lim')))
}