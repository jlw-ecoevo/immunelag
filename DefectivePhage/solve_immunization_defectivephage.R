

#Early infection dynamics

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
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- (v*R/(z+R) - delta*V - nu*delta*Vdefective - w)*U
    dI <- delta*U*V - gam*I - w*I 
    dD <- (v*R/(z+R) - w - delta*V)*D + nu*delta*Vdefective*U + phi*Id
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0*(1-frac_defective)-V) + beta*(1-frac_defective)*gam*I - delta*(U+D)*V
    dVdefective <- w*(v0*(frac_defective)-Vdefective) + beta*(frac_defective)*gam*I - delta*(U+D)*Vdefective
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV,dVdefective))
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
    dU <- (v*R/(z+R) - delta*V - nu*delta*Vdefective - w)*U
    dI <- delta*U*V - gam*I - w*I 
    dD <- (v*R/(z+R) - w)*D + nu*delta*Vdefective*U
    dId <-  0
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0*(1-frac_defective)-V) + beta*(1-frac_defective)*gam*I - delta*(U+D)*V
    dVdefective <- w*(v0*(frac_defective)-Vdefective) + beta*(frac_defective)*gam*I - delta*(U+D)*Vdefective
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV,dVdefective))
  }) # end with(as.list ...
}



getParameters <- function(w=0.3,
                          r0=350,
                          e=5e-7,
                          v=2,
                          z=1,
                          delta=1e-7,
                          mu=1e-7,
                          gam=4/3,
                          beta=80,
                          v0=0,
                          phi=4/3,
                          xi=1/6,
                          kappa=0.05,
                          frac_defective=0.1,
                          nu=mu/frac_defective){
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
                  kappa=kappa,
                  frac_defective=frac_defective,
                  nu=nu)
  return(parameters)
}

# Helper function to change initial conditions
getInitial <- function(R=350,
                       U=1e8,
                       I=0,
                       D=100,
                       Id=0,
                       SM=100,
                       V=0,
                       Vdefective=0){
  state <- c(R=R,
             U=U,
             I=I,
             D=D,
             Id=Id,
             SM=SM,
             V=V,
             Vdefective=Vdefective)
  return(state)
}

plotOut <- function(out,main=""){
  par(mar=c(5.1, 5, 4.1, 1))
  plot(out[,"time"],out[,"U"],type="l",
       log="y",ylim=c(1,1e11),lwd=2,main=main,
       xlab="Time (hours)",ylab="Density",
       cex.axis=1.5,cex.lab=1.5,cex.main=1)
  lines(out[,"time"],out[,"I"],col="black",lty=2,lwd=3)
  lines(out[,"time"],out[,"D"],col=pal_npg("nrc")(6)[4],lwd=5)
  lines(out[,"time"],out[,"Id"],col=pal_npg("nrc")(6)[4],lty=2,lwd=3)
  lines(out[,"time"],out[,"SM"],col=pal_npg("nrc")(6)[3],lwd=5)
  lines(out[,"time"],out[,"V"],col=pal_npg("nrc")(6)[1],lwd=5)
  par(mar=c(5.1, 4.1, 4.1, 2.1))
}



times <- seq(0, 24, by = 1)

parameters <- getParameters(kappa=0.01,phi=10)
state <- getInitial(V=100)
out <- ode(y = state, 
           times = times, 
           func = autoimmunitySystem, 
           parms = parameters,
           method="lsoda")

parameters <- getParameters(kappa=0.01,phi=10)
state <- getInitial(V=100)
out_nl <- ode(y = state, 
           times = times, 
           func = autoimmunitySystemNoLag, 
           parms = parameters,
           method="lsoda")

parameters <- getParameters(kappa=0.01,phi=10)
state <- getInitial(V=100,D=0,SM=1)
out_im <- ode(y = state, 
           times = times, 
           func = autoimmunitySystem, 
           parms = parameters,
           method="lsoda")

parameters <- getParameters(kappa=0.01,phi=10)
state <- getInitial(V=100,D=0,SM=1)
out_nl_im <- ode(y = state, 
              times = times, 
              func = autoimmunitySystemNoLag, 
              parms = parameters,
              method="lsoda")

setwd("~/immunelag/DefectivePhage/")
pdf(paste0("ImmunizationLag_defectivephage.pdf"),width=10,height=10)
par(mfrow=c(2,2))
plotOut(out,main=expression("Immunization With Lag"))
plotOut(out_nl,main=expression("Immunization Without Lag"))
plotOut(out_im,main=expression("Immunization With Lag (Non-Immune Start)"))
plotOut(out_nl_im,main=expression("Immunization Without Lag (Non-Immune Start)"))
par(mfrow=c(1,1))
dev.off()

