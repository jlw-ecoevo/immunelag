
# JLW -2020

# Early infection dynamics during an outbreak (24 hrs w/ susceptible initial population)
# Panels A, B, C

# Load Packages
library(deSolve)
library(ggsci)
library(scales)

# Lag System
autoimmunitySystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- (v*R/(z+R) - delta*V - w)*U
    dI <- (1-mu)*delta*U*V - gam*I - w*I
    dD <- (v*R/(z+R) - w - delta*V)*D + mu*delta*V*U + phi*Id
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV))
  }) # end with(as.list ...
}

# System of ODEs w/out lag
autoimmunitySystemNoLag <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- (v*R/(z+R) - delta*V - w)*U
    dI <- (1-mu)*delta*U*V - gam*I - w*I
    dD <- (v*R/(z+R) - w)*D + mu*delta*V*U
    dId <-  0
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV))
  }) # end with(as.list ...
}

# Helper function to change parameters
getParameters <- function(w=0.3,
                r0=350,
                e=5e-7,
                v=2,
                z=1,
                mu=1e-7,
                delta=1e-7,
                gam=4/3,
                beta=80,
                v0=0,
                phi=4/3,
                xi=1/6,
                kappa=0.05){
  parameters <- c(w=w,
                  r0=r0,
                  e=e,
                  v=v,
                  z=z,
                  delta=delta,
                  mu=mu,
                  gam=gam,
                  beta=beta,
                  v0=v0,
                  phi=phi,
                  xi=xi,
                  kappa=kappa)
  return(parameters)
}

# Helper function to change initial conditions
getInitial <- function(R=350,
                       U=1e8,
                       I=0,
                       D=100,
                       Id=0,
                       SM=100,
                       V=0){
  state <- c(R=R,
             U=U,
             I=I,
             D=D,
             Id=Id,
             SM=SM,
             V=V)
  return(state)
}

# Helper function to plot model output
plotOut <- function(out,main=""){
  par(mar=c(5.1, 5, 4.1, 1))
  plot(out[,"time"],out[,"U"],type="l",
       log="y",ylim=c(1,1e11),lwd=2,main=main,
       xlab="Time (hours)",ylab="Density",
       cex.axis=1.5,cex.lab=1.5,cex.main=2)
  lines(out[,"time"],out[,"I"],col="black",lty=2,lwd=3)
  lines(out[,"time"],out[,"D"],col=pal_npg("nrc")(6)[4],lwd=5)
  lines(out[,"time"],out[,"Id"],col=pal_npg("nrc")(6)[4],lty=2,lwd=3)
  lines(out[,"time"],out[,"SM"],col=pal_npg("nrc")(6)[3],lwd=5)
  lines(out[,"time"],out[,"V"],col=pal_npg("nrc")(6)[1],lwd=5)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}

# Helper function to plot model output w/ CRISPR defended combined
plotOutCombine <- function(out,main=""){
  par(mar=c(5.1, 5, 4.1, 1))
  plot(out[,"time"],out[,"U"],type="l",
       log="y",ylim=c(1,1e11),lwd=2,main=main,
       xlab="Time (hours)",ylab="Density",
       cex.axis=1.5,cex.lab=1.5,cex.main=2)
  lines(out[,"time"],out[,"I"],col="black",lty=2,lwd=3)
  lines(out[,"time"],out[,"D"]+out[,"Id"],col=pal_npg("nrc")(6)[4],lwd=5)
  lines(out[,"time"],out[,"SM"],col=pal_npg("nrc")(6)[3],lwd=5)
  lines(out[,"time"],out[,"V"],col=pal_npg("nrc")(6)[1],lwd=5)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}


# Time range and initial conditions for numerical soln
times <- seq(0, 24, by = 1)
state <- getInitial(V=100)

# Run model w/ no, short, long lags
parameters <- getParameters(v0=100,kappa=0.01,phi=10)
out_immunization_longlag <- ode(y = state, 
           times = times, 
           func = autoimmunitySystem, 
           parms = parameters,
           method="lsoda")
parameters <- getParameters(v0=100,kappa=0.01,phi=1e3)
out_immunization_shortlag <- ode(y = state, 
           times = times, 
           func = autoimmunitySystem, 
           parms = parameters,
           method="lsoda")
parameters <- getParameters(v0=100,kappa=0.01,phi=0)
out_immunization_nolag <- ode(y = state, 
           times = times, 
           func = autoimmunitySystemNoLag, 
           parms = parameters,
           method="lsoda")

#save model output
setwd("~/immunelag/Fig4")
save(out_immunization_nolag,
     out_immunization_shortlag,
     out_immunization_longlag,
     file = "immunizationlag.RData")

# Plot
setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()

# Plot
setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLagCombine.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOutCombine(out_immunization_nolag,main=expression("Outbreak, No Lag"))
plotOutCombine(out_immunization_shortlag,main=expression("Outbreak, Short Lag"))
plotOutCombine(out_immunization_longlag,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()


# Plot
setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_slides.pdf"),width=4,height=8)
par(mfrow=c(2,1))
plotOut(out_immunization_nolag,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_longlag,main=expression("Outbreak, Lag"))
par(mfrow=c(1,1))
dev.off()


# Very costly surface mutant (kappa = 0.1, otherwise same as above)

parameters <- getParameters(v0=100,kappa=0.1,phi=10)
out_immunization_longlag_costlysm <- ode(y = state, 
                                times = times, 
                                func = autoimmunitySystem, 
                                parms = parameters,
                                method="lsoda")
parameters <- getParameters(v0=100,kappa=0.1,phi=1e3)
out_immunization_shortlag_costlysm <- ode(y = state, 
                                 times = times, 
                                 func = autoimmunitySystem, 
                                 parms = parameters,
                                 method="lsoda")
parameters <- getParameters(v0=100,kappa=0.1,phi=0)
out_immunization_nolag_costlysm <- ode(y = state, 
                              times = times, 
                              func = autoimmunitySystemNoLag, 
                              parms = parameters,
                              method="lsoda")


setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_costly.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag_costlysm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlysm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlysm,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()

parameters <- getParameters(v0=100,kappa=0.2,phi=10)
out_immunization_longlag_costlyx2sm <- ode(y = state, 
                                             times = times, 
                                             func = autoimmunitySystem, 
                                             parms = parameters,
                                             method="lsoda")
parameters <- getParameters(v0=100,kappa=0.2,phi=1e3)
out_immunization_shortlag_costlyx2sm <- ode(y = state, 
                                              times = times, 
                                              func = autoimmunitySystem, 
                                              parms = parameters,
                                              method="lsoda")
parameters <- getParameters(v0=100,kappa=0.2,phi=0)
out_immunization_nolag_costlyx2sm <- ode(y = state, 
                                           times = times, 
                                           func = autoimmunitySystemNoLag, 
                                           parms = parameters,
                                           method="lsoda")

setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_costlyx2.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag_costlyx2sm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlyx2sm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlyx2sm,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()



parameters <- getParameters(v0=100,kappa=0.4,phi=10)
out_immunization_longlag_costlyx4sm <- ode(y = state, 
                                             times = times, 
                                             func = autoimmunitySystem, 
                                             parms = parameters,
                                             method="lsoda")
parameters <- getParameters(v0=100,kappa=0.4,phi=1e3)
out_immunization_shortlag_costlyx4sm <- ode(y = state, 
                                              times = times, 
                                              func = autoimmunitySystem, 
                                              parms = parameters,
                                              method="lsoda")
parameters <- getParameters(v0=100,kappa=0.4,phi=0)
out_immunization_nolag_costlyx4sm <- ode(y = state, 
                                           times = times, 
                                           func = autoimmunitySystemNoLag, 
                                           parms = parameters,
                                           method="lsoda")

setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_costlyx4.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag_costlyx4sm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlyx4sm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlyx4sm,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()


parameters <- getParameters(v0=100,kappa=0,phi=10)
out_immunization_longlag_nocostsm <- ode(y = state, 
                                           times = times, 
                                           func = autoimmunitySystem, 
                                           parms = parameters,
                                           method="lsoda")
parameters <- getParameters(v0=100,kappa=0,phi=1e3)
out_immunization_shortlag_nocostsm <- ode(y = state, 
                                            times = times, 
                                            func = autoimmunitySystem, 
                                            parms = parameters,
                                            method="lsoda")
parameters <- getParameters(v0=100,kappa=0,phi=0)
out_immunization_nolag_nocostsm <- ode(y = state, 
                                         times = times, 
                                         func = autoimmunitySystemNoLag, 
                                         parms = parameters,
                                         method="lsoda")

setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_nocost.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag_nocostsm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_nocostsm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_nocostsm,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()



setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag_comparecost.pdf"),width=12,height=16)
par(mfrow=c(4,3))
plotOut(out_immunization_nolag_nocostsm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_nocostsm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_nocostsm,main=expression("Outbreak, Long Lag"))

plotOut(out_immunization_nolag_costlysm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlysm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlysm,main=expression("Outbreak, Long Lag"))

plotOut(out_immunization_nolag_costlyx2sm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlyx2sm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlyx2sm,main=expression("Outbreak, Long Lag"))

plotOut(out_immunization_nolag_costlyx4sm,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag_costlyx4sm,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag_costlyx4sm,main=expression("Outbreak, Long Lag"))
par(mfrow=c(4,3))
dev.off()


