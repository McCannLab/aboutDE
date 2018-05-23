#' Simple numerical integration
#'
#' MyRK4 function performs a numerical integration of function FUN
#' using the Runge-Kutta at order 4 methods. (no optimization!)
#' a R function describing the set of ODE and a initial condition must be provided
#'
#'
#' @author
#' Kevin Cazelles

#' @param Xi: vector containing initial conditions.
#' @param ti: initial time, 0 by default.
#' @param tf: final time, 1 by dafult.
#' @param d_t: time increment, 0.001 by default.
#' @param ...: additionnal argument to be passed to FUN words. Not \code{NA_character_}.
#'
#' @export
#' @examples
# systLin <- function(X, beta){
#     Y <- matrix(0,ncol=2)
#     Y[1] <- beta[1,1]*X[1]+beta[1,2]*X[2]
#     Y[2] <- beta[2,1]*X[1]+beta[2,2]*X[2]
#     return(Y)
# }
# A <- matrix(c(0,1,-1,0),2)
# res <-RK4(ti=0, tf=10, d_t=0.01, FUN=systLin, Xi=c(0.8,0.8), beta=A)
# plot(res)

RK4 <- function(FUN, Xi, ti=0, tf=1, d_t=0.001, ...){
    ## Time sequence
    seqt <- seq(ti,tf,d_t)
    nt <- length(seqt)
    ## Values
    X <- matrix(0, ncol=length(Xi), nrow=nt)
    X[1,] <- Xi
    ## Loop
    for (i in 1:(nt-1)){
        k1 <- FUN(X[i,],...)
        k2 <- FUN(X[i,]+.5*d_t*k1,...)
        k3 <- FUN(X[i,]+.5*d_t*k2,...)
        k4 <- FUN(X[i,]+d_t*k3,...)
        X[i+1,] <- X[i,]+ d_t*(1/6)*(k1+2*k2+2*k3+k4)
    }

    data.frame(seqt,X)
}
