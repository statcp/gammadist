#' Construct lower tolerance limit, upper tolerance limit and tolerance interval for a
#' proportion of future measurements from the gamma distribution.
#'
#' @param x observations of a gamma distribution.
#' @param alpha value of alpha, such that a (beta, 1-alpha) lower tolerance limit, a (beta, 1-alpha)
#' upper tolerance limit or a (beta, 1-alpha) tolerance interval is returned. If alpha is not specified, the default
#' value is 0.05.
#' @param gamma value of gamma (the proportion of future measurements), such that a (gamma, 1-alpha) lower
#' tolerance limit, a (gamma, 1-alpha) upper tolerance limit or a (gamma, 1-alpha) tolerance interval
#' is returned. If gamma is not specified, the default value is 0.95.
#' @param B number of realizations of the GPQs. If B is not specified, the default value is 2000.
#'
#' @details Assume that X_1, ..., X_n from a gamma distribution with culmulative distribution function F are available. 
#' Consider two statistics TL = TL(X_1, ..., X_n) and TU = PU(X_1, ..., X_n). 
#' Let gamma and alpha be two constants and 0 < gamma, alpha < 1. If TL and TU
#' are determined such that
#' P(F(TU) - F(TL) >= gamma) = 1-alpha,
#' then (TL,TU) is called a (gamma, 1-alpha) two-sided tolerance interval. When TU equals infinity, TL is called
#' a (gamma, 1-alpha) lower tolerance limit. 
#' Similarly, when TL equals -infinity, TU is called a (gamma, 1-alpha) upper tolerance limit. 
#'
#' @return  tollimits returns a dataframe containing one vector where the first two elements make for the (gamma, 1-alpha) tolerance
#'  interval, the third element is the 
#' (gamma, 1-alpha) lower tolerance limit and the last is the (gamma, 1-alpha) 
#'  upper tolerance limit.
#' 
#' @export
#'
#' @references Chen, P. and Ye, Z., 2017. "Approximate Statistical Limits For A Gamma Distribution."
#' Journal of Quality Technology, 49(1), 64-77.
#'
#'Wang, B. X. and Wu, F. (2018), “Inference on the Gamma Distribution,” 
#'Technometrics, 60(2), 235–244.
#'
#' Krishnamoorthy, K. Mathew, T. Mukherjee, S. 2008. "Normal-Based Methods for a Gamma Distribution:
#' Prediction and Tolerance Intervals and Stress-Strength Reliability". 
#' Technometrics 50(1), 69-78.
#' 
#' @examples x = rGamma(100, shape = 3, rate = 1)
#' tollimits(x, alpha = 0.05, gamma = 0.99, B=5000)


tollimits <- function(x, alpha = 0.05, gamma = 0.99, B=5000){
  n <- length(x)
  w <- x^(1/3)
  wbar <- mean(w)
  sw <- sd(w)
  c <- sqrt((n-1)*qchisq(gamma,df=1,ncp=1/n)/qchisq(alpha,n-1))
  gpq.all <- pargpq (x,B)
  tu.all <- mapply(function(shape,scale) qgamma(gamma,shape=shape,scale=scale), shape=gpq.all$shape,scale=gpq.all$scale)
  tl.all <- mapply(function(shape,scale) qgamma(1-gamma,shape=shape,scale=scale), shape=gpq.all$shape,scale=gpq.all$scale)
  tol <- c((wbar-c*sw)^3,(wbar+c*sw)^3,quantile(tl.all,alpha),quantile(tu.all,1-alpha))
  return(data.frame(tol,row.names=c('low-int','up-int','low-lim','up-lim')))
}

