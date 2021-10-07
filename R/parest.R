#' Estimate the parameters of the gamma distribution.
#'
#' @param x gamma distributed observations.
#'
#' @return parest gives the estimates of the shape parameter, the rate parameter and and the scale parameter.
#'
#' @export
#'
#' @references Ye, Z.-S. and Chen, N. (2017), “Closed-form estimators for the gamma distribution derived from likelihood equations,” The American Statistician, 71(2), 177–181.
#'
#'Louzada, F., Ramos, P. L., and Ramos, E. (2019), “A note on bias of closed-form estimators for the gamma distribution derived from likelihood equations,” The American Statistician, 73(2), 195–199.
#'
#'
#' @examples x <- rGamma(1000, shape = 2, rate = 1)
#' parest(x)
#'
#' x <- c(1.6, 2.7, 1.4, 2.4, 1.3, 1.8, 2.3, 1.1)
#' parest(x)


parest <- function(x){
  n <- length(x)
  shape0 <- (n*sum(x)) / (n*sum(x*log(x)) - sum(log(x))*sum(x))
  rate0 <- n^2 / (n*sum(x*log(x)) - sum(log(x))*sum(x))
  scale0 <- 1/rate0
  shape <- (n-1)*shape0/(n+2)
  rate <- (n-1)*rate0/(n+2)
  scale <- n*scale0/(n-1)
  parameters <- list("shape" =  shape, "rate" = rate, "scale" = scale)
  return(parameters)
}
