#' Generate the realizations of the GPQs for the parameters of the gamma distribution.
#'
#' @param x observations of a gamma distribution.
#' @param B number of realizations of the GPQs. If B is not specified, the default value is 2000.
#'
#' @details Assume that the data x_1, ..., x_n from a gamma distribution are available. 
#' The GPQs for the gamma parameters have been developed in Wang and Wu (2018). 
#' Because the exact distributions of the GPQs are difficult to derive, 
#' this pargpq function gives realizations of them using the Monte Carlo procedures.

#' @return pargpq returns a dataframe which contains a list of the variables shape, scale and rate. Each variable has B realizations of the corresponding GPQ.
#' @export
#'
#' @references Wang, B. X. and Wu, F. (2018), “Inference on the Gamma Distribution,” Technometrics, 60(2), 235–244.

#' @examples x <- rGamma(100, shape = 3, rate = 1)
#' pargpq(x)
#' pargpq(x, B=10000)


pargpq <- function(x,B=2000){
  n <- length(x)
  xtilde <- (prod(x))^(1/n); xbar <- mean(x);
  t0 <- log(xtilde/xbar)
  cul <- function(k){
    i <- 2:5
    c(psigamma(k)+log(n)-psigamma(n*k),psigamma(k,i-1)/n^(i-1)-psigamma(n*k,i-1))
  }
  t <- function(k,p){
    cc <- cul(k)
    ccd <- cc/cc[2]^(1:5/2)
    zp <- qnorm(p)
    z <- zp+ccd[3]*(zp^2-1)/6+ccd[4]*(zp^3-3*zp)/24-ccd[3]^2*(2*zp^3-5*zp)/36+ccd[5]*(zp^4-6*zp^2+3)/120-ccd[3]*ccd[4]*(zp^4-5*zp^2+2)/24+ccd[3]^3*(12*zp^4-53*zp^2+17)/324
    cc[1]+cc[2]^0.5*z-t0
  }
  solvek <- function(p){
    tryCatch(uniroot(t,p=p,interval=c(10^-2,100))$root,error=function(e) NA)
  }
  shape <- sapply(runif(B),solvek,simplify=T)
  scale <- 2*n*xbar/rchisq(B,2*n*shape)
  rate <- 1/scale
  return(data.frame(shape,scale,rate))
}



