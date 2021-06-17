
# JLW - 2020

# Analysis of experimental data and comparison with model expectations

# Load packages
library(deSolve)
library(ggsci)
library(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)

parameters <- c(w=0.3,
                r0=350,
                e=5e-7,
                v=2,
                z=1,
                delta=1e-7,
                gam=4/3,
                beta=80,
                v0=0,
                phi=0,
                k=0.1)

state <- c(
  R=500,
  V=0,
  C=1e6,
  L=0,
  M=1e6
)

times <- seq(0, 24, by = 1)

autoimmunitySystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <-  - e*v*R*(C+M)/(z+R)
    dC <- (v*R/(z+R) - delta*V)*C + phi*L
    dL <-  delta*C*V - phi*L
    dM <- ((1-k)*v*R/(z+R))*M
    dV <- - delta*C*V
    
    # return the rate of change
    list(c(dR,dV,dC,dL,dM))
  }) # end with(as.list ...
}

plotOut <- function(out,main=""){
  par(mar=c(5.1, 5, 4.1, 1))
  plot(out[,"time"],out[,"C"],col=pal_npg("nrc")(6)[4],
       lwd=5,main=main,type="l",log="y",
       xlab="Time (hours)",ylab="Density",
       cex.axis=1.5,cex.lab=1.5,cex.main=2,
       ylim=c(1,1e9))
  lines(out[,"time"],out[,"M"],col=pal_npg("nrc")(6)[3],lwd=5)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

out <- ode(y = state, 
           times = times, 
           func = autoimmunitySystem, 
           parms = parameters,
           method="lsoda")

plotOut(out)

Yf <- tail(out[,"C"],n=1)
Xf <- tail(out[,"M"],n=1)
Xprop <- Xf/(Xf+Yf)
rel_fit <- (Xprop*0.5)/((1-Xprop)*0.5)
rel_fit
