#==============================================================================
#==============================================================================
#    2-dimensional systems
#
# This script exemplifies how to simulate 2 dimensional ordinary 
# differential equations using R.
#
# code written by: Torbjörn Säterberg 180523 
#
#==============================================================================
#==============================================================================

# load packages
library(deSolve)

# clear environment
rm(list=ls())

#===============================================================================
# Model 1
# 
# Standard Lotka-Volterra model
#
# dRdt = r * R  - alpha * R * C
# dCdt = e * alpha * R * C - m * C 
#===============================================================================

#-------------------#
# 1) model function #
#-------------------#

LVmod1 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dR    <- r * R - alpha * R * C
    dC    <- alpha * R * C * e - m * C 
    
    return(list(c(dR, dC)))
  })
}

#-------------------#
# 2) Zero isoclines #
#-------------------#
# dRdt=0 =>  C=r/alpha
# dCdt=0 =>  R=m/e*alpha

LVmod1Zeroisocline <- function(pars) {
  with(as.list(pars),{
    
    # Resource isocline
    plot(c(0,2*m/(alpha*e)),c(r/alpha,r/alpha),col="green",type="l",xlab="Resource density",ylab="Consumer density")
    
    # Consumer isocline
    lines(c(m/(e*alpha),m/(e*alpha)),c(0,2*r/alpha),col="red",main="Zero isoclines")  
    legend("topright", c("dRdt=0", "dCdt=0"), col =c("green","red"), lty = 1)
  })
}

#------------------------------------------------#
# 3) Set parameters and simulate data            #
#------------------------------------------------#

pars    <- c(alpha = 0.2,    # attack rate
             r = 1.0,    # per capita growth rate of prey
             m = 0.2 ,   # per capita mortality rate of consumer
             e = 0.5,    # conversion efficiency
             K  = 10)     # carrying capacity

yini    <- c(R = 1, C = 2)
times   <- seq(0, 200, by = 0.01)

# simulate model
out     <- lsoda(func = LVmod1, y = yini, parms = pars, times = times)

#----------------#
# 4) Plot model  #
#----------------#

# Subplot set-up
par(mfrow=c(1,3))
par(mar=c(c(4, 4, 2, 1)))

# Time trajectories
matplot(out[,"time"], out[,2:3], type = "l", xlab = "time", ylab = "Conc",
        main = "Lotka-Volterra", lwd = 2)
legend("topright", c("Resource", "Consumer"), col =1:2, lty = 1:2)

# Phase space plot
matplot(out[,"R"],out[,"C"],type="l",xlab="Resource density",
        ylab="Consumer density", main="Phase space plot")

# Zero isoclines
LVmod1Zeroisocline(pars)


#===============================================================================
# Model 2
#
# Lotka-volterra with logistic growth of resource.
#
# dRdt = r * R * (1 - R/K) - alpha * R * C
# dCdt = alpha * R * C * e - m * C 
#===============================================================================

#-------------------#
# 1) model function #
#-------------------#

LVmod2 <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    
    dR    <- r * R * (1 - R/K) - alpha * R * C
    dC    <- e * alpha * R * C  - m * C 
    
    return(list(c(dR, dC)))
  })
}

#-------------------# 
# 2) Zero isoclines #
#-------------------#
# dRdt=0 =>  C=(r/alpha)*(1-R/K)
# dCdt=0 =>  R=m/(e*alpha)

LVmod2Zeroisocline <- function(pars) {
  with(as.list(pars),{
    
    # Resource isocline
    plot(c(0,K),c(r/alpha,0),col="green",type="l",
         xlab="Resource density",ylab="Consumer density",
         main="Zero isoclines")
    
    # Consumer isocline
    lines(c(m/(e*alpha),m/(e*alpha)),c(0,2*r/alpha),col="red")  
    legend("topright", c("dRdt=0", "dCdt=0"), col =c("green","red"), lty = 1)
  })
}

#------------------------------------------------#
# 3) Set parameters and simulate data            #
#------------------------------------------------#

pars    <- c(alpha = 0.2,    # attack rate
             r = 1.0,    # per capita growth rate of prey
             m = 0.2 ,   # per capita mortality rate of consumer
             e = 0.5,    # conversion efficiency
             K  = 10)     # carrying capacity

yini    <- c(R = 1, C = 2)
times   <- seq(0, 100, by = 0.01)

# simulate model
out     <- lsoda(func = LVmod2, y = yini, parms = pars, times = times)

#----------------#
# 4) Plot model  #
#----------------#

# Subplot set-up
par(mfrow=c(1,3))
par(mar=c(c(4, 4, 2, 1)))

# Time trajectories
matplot(out[,"time"], out[,2:3], type = "l", xlab = "Time", ylab = "Density",
        main = "Time trajectories", lwd = 2)
legend("topright", c("Resource", "Consumer"), col =1:2, lty = 1:2)

# Phase space plot
matplot(out[,"R"],out[,"C"],type="l",xlab="Resource density",
        ylab="Consumer density",main="Phase space plot")

# Zero isoclines
LVmod2Zeroisocline(pars)

