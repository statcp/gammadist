#' Generate gamma distributed random variables.
#'
#' @param n number of observations.
#' @param shape value of the shape parameter (> 0).
#' @param rate value of the rate parameter (> 0).
#' @param scale value of the scale parameter (> 0).
#'
#' @details If rate and scale are ommitted, it assumes the default value of 1.
#'
#' @return rGamma(n, shape, rate, scale) generates n gamma distributed random numbers with parameters
#' shape and rate.
#' @export
#'
#' @references Liu, C., Martin, R., and Syring, N. (2017), “Efficient simulation from a gamma distribution with
#' small shape parameter,” Computational Statistics, 32(4), 1767–1775.
#'
#' @examples rGamma(100, shape = 0.5, rate = 3)
#'
#' rGamma(1000, shape = 2, scale = 2)


rGamma <- function(n, shape, rate = 1, scale = 1/rate){
  rate = 1/scale
  X = numeric(n)
  if (shape >= 1) {
    X = rgamma(n, shape = shape, rate = rate)
    return(X)
  } else{
    w <- shape/2.71828182845905/(1-shape); lambda <- 1/shape-1; r <- 1/(1+w)
    eta <- function(z) ifelse(z >= 0, exp(-z), w *lambda*exp(lambda*z))
    h <- function(z) exp(-z-exp(-z/shape))
    for (i in 1:n){
      repeat{
        z <- ifelse(runif(1) <=r,rexp(1),-rexp(1,lambda))
        if (h(z)/eta(z)>runif(1)){
          X[i] <- exp(-z/shape)/rate
          break
        }
      }
    }
    return(X)
  }
}
