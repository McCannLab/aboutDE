library(odeintr)
library(graphicutils)

# linear system with ODEint



intLinear <- function(params, init = c(1,1), mxt=10, dt = 0.01, beta = beta1, seqx = seq(-2,2,.41), seqy = seq(-2,2,.41)) {

  systLin <- function(X, t){
     Y <- matrix(0,ncol=2)
     Y[1] <- beta[1,1]*X[1]+beta[1,2]*X[2]
     Y[2] <- beta[2,1]*X[1]+beta[2,2]*X[2]
     return(Y)
  }

  #
  obs <- function(x, t) c(sp1 = x[1], sp2 = x[2])
  x <- integrate_sys(systLin, init, mxt, dt, observer = obs)

  ##
  par(mar=c(1,1,1,1))
  plot0(c(-2,2),c(-2,2))
   vecfield2d(coords=expand.grid(seqx, seqy), FUN=systLin,
      cex.x=0.35, cex.arr=0.25,
      border=NA,cex.hh=1, cex.shr=0.6, col=8, add = T)
    points(0, 0, pch =19, cex = 2)
    points(init[1], init[2], col = "#207966", pch =19, cex = 2)
    lines(x$sp1, x$sp2, col = "#207966", lwd = 2)
  abline(v=0,h=0)

}

beta1 <- matrix(c(-1, 0, 0, -1), 2)
beta1b <- matrix(c(-1, 0, 0, 1), 2)
beta1c <- matrix(c(1, 0, 0, 1), 2)
##
beta2 <- matrix(c(0,-1,1,0), 2)
beta2b <- matrix(c(0,1,1,0), 2)
beta2c <- matrix(c(0,-1,-1,0), 2)
##
beta3 <- matrix(c(0.1,-1,1,0.1), 2)
beta4 <- matrix(c(-0.1,-1,1,-0.1), 2)

intLinear()
intLinear(beta = beta1b)
intLinear(beta = beta1c)
intLinear(beta = beta2)
intLinear(beta = beta2b)
intLinear(beta = beta2c)

intLinear(beta = beta3)
intLinear(beta = beta4, mxt =100)





# Solve sys with RK4.R
source("R/RK4.R")
systLin <- function(X, beta, t){
   Y <- matrix(0,ncol=2)
   Y[1] <- beta[1,1]*X[1]+beta[1,2]*X[2]
   Y[2] <- beta[2,1]*X[1]+beta[2,2]*X[2]
   return(Y)
}

ti = 0
tf = 10
d_t = 0.01
out <- RK4(ti=ti, tf=tf, d_t=d_t, FUN=systLin, Xi=c(0.8,0.8), beta=beta3)

vecfield2d(coords=expand.grid(seqx, seqy), FUN=systLin,
   args=list(beta=beta3), cex.x=0.35, cex.arr=0.25,
   border=NA,cex.hh=1, cex.shr=0.6, col=8)
   lines(out$X1, out$X2)
# abline(a=c(0,0), b=c(.1, 10))
abline(v=0,h=0)
