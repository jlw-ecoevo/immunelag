

#Early infection dynamics during an outbreak (24 hrs w/ susceptible initial population)

library(deSolve)
library(ggsci)
library(scales)

autoimmunitySystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Iv <- max(Iv,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- ((1-2*mu/alpha)*v*R/(z+R) - delta*V - w)*U
    dI <- delta*U*V - gam*I - w*I - mu*Iv
    dIv <- delta*U*V + gam*log(beta)*Iv - (gam*exp(xi*log(beta))/xi)*I - mu*(Iv^2)/(I+1e-16)
    dD <- ((1-2*mu/alpha)*v*R/(z+R) - w - delta*V)*D + mu*Iv + phi*Id
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dIv,dD,dId,dSM,dV))
  }) # end with(as.list ...
}

autoimmunitySystemNoLag <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Iv <- max(Iv,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- ((1-2*mu/alpha)*v*R/(z+R) - delta*V - w)*U
    dI <- delta*U*V - gam*I - w*I - mu*Iv
    dIv <- delta*U*V + gam*log(beta)*Iv - (gam*exp(xi*log(beta))/xi)*I - mu*(Iv^2)/(I+1e-16)
    dD <- ((1-2*mu/alpha)*v*R/(z+R) - w)*D + mu*Iv 
    dId <-  0
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dIv,dD,dId,dSM,dV))
  }) # end with(as.list ...
}

getParameters <- function(w=0.3,
                r0=350,
                e=5e-7,
                v=2,
                z=1,
                delta=1e-7,
                mu=1e-6,
                alpha=0.4,
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
                  alpha=alpha,
                  gam=gam,
                  beta=beta,
                  v0=v0,
                  phi=phi,
                  xi=xi,
                  kappa=kappa)
  return(parameters)
}

getInitial <- function(R=350,
                       U=1e8,
                       I=0,
                       Iv=0,
                       D=100,
                       Id=0,
                       SM=100,
                       V=0){
  state <- c(R=R,
             U=U,
             I=I,
             Iv=Iv,
             D=D,
             Id=Id,
             SM=SM,
             V=V)
  return(state)
}

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


times <- seq(0, 24, by = 1)
state <- getInitial(V=100)

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

setwd("~/immunelag/Fig4")
save(out_immunization_nolag,
     out_immunization_shortlag,
     out_immunization_longlag,
     file = "immunizationlag.RData")

setwd("~/immunelag/Fig4")
pdf(paste0("ImmunizationLag.pdf"),width=12,height=4)
par(mfrow=c(1,3))
plotOut(out_immunization_nolag,main=expression("Outbreak, No Lag"))
plotOut(out_immunization_shortlag,main=expression("Outbreak, Short Lag"))
plotOut(out_immunization_longlag,main=expression("Outbreak, Long Lag"))
par(mfrow=c(1,3))
dev.off()


# Very costly surface mutant

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

