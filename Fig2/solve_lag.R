
# JLW - 2020

# Panel B for Figure 2 (numerical solutions for lag model)

# Load packages
library(deSolve)
library(wesanderson)

# System of differential equations
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
    dIv <- delta*U*V - w*Iv - mu*((Iv^2)/(I+1e-16))
    dD <- ((1-2*mu/alpha)*v*R/(z+R) - w - delta*V)*D + mu*Iv + phi*Id
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dIv,dD,dId,dSM,dV))
  }) # end with(as.list ...
}


# Given a pair of parameters phi and v0, returns state of system at 1e7 hours (presumed equilibrium)
phi.vs.v0 <- function(phi,v0,state,parameters){
  parameters["v0"] <- v0
  parameters["phi"] <- phi
  times <- seq(0, 1e7, by = 1e6)
  out <- ode(y = state, 
             times = times, 
             func = autoimmunitySystem, 
             parms = parameters,
             method="lsoda")
  return(list(D = (tail(out[,"D"],n=1) < 1),
              SM = (tail(out[,"SM"],n=1) > 1e8)))
}

# Runs analysis over several SM costs
sweepKappa <- function(k){
  parameters <- c(w=0.3,
                  r0=350,
                  e=5e-7,
                  v=2,
                  z=1,
                  delta=1e-7,
                  mu=1e-6,
                  alpha=0.4,
                  gam=4/3,
                  beta=80,
                  v0=1e5,
                  phi=4/3,
                  kappa=k)
  
  state <- c(R=350,
             U=0,
             I=0,
             Iv=0,
             D=1e8,
             Id=0,
             SM=100,
             V=0)
  
  #phi_range <-  rev(seq(0.1,5,0.25))
  phi_range <- (10^seq(-5,3.5,0.1))
  v0_range <- 10^seq(0,12,0.1)
  pv_combos <- expand.grid(phi_range,v0_range)
  names(pv_combos) <- c("phi","v0")
  x <- mapply(phi.vs.v0,phi=pv_combos$phi,v0=pv_combos$v0,MoreArgs = list(state=state,parameters=parameters))
  
  pv_combos$Invasion <- (x["D",]==T & x["SM",]==T)*1
  i_mat <- reshape(pv_combos, idvar = "phi", timevar = "v0", direction = "wide")
  
  return(i_mat[,-1])
}

# Numeric solutions for model varying kappa, phi, and v0 (takes a long time to run)
kappas <- c(0.001,0.01,0.1)
x <- lapply(kappas,sweepKappa)
setwd("~/immunelag/Fig2")
save(x,file="kappa_sweep.RData")

setwd("~/immunelag/Fig2")
load("kappa_sweep.RData")

# Plot is all in a nice contour plot
phi_range <- (10^seq(-5,3.5,0.1))
v0_range <- 10^seq(0,12,.1)
setwd("~/immunelag/Fig2")
pdf(paste0("InvasionLag_contour.pdf"),width=8,height=5)
par(mar=c(5.1, 5, 1, 0.2))
filled.contour(x=log10(phi_range),
               y=log10(v0_range),
               z=(as.matrix(x[[2]])),
               xlab=expression("Recovery Rate from Immune Lag (" * log[10](phi) * ")"),
               ylab=expression("Envrionmental Viral Pool (" * log[10](v0) * ")"),
               levels=c(0,0.5,1),
               col = wes_palette("Moonrise3",2),
               cex.lab=1.5,
               cex.axis=1.5,
               plot.axes = { contour(x=log10(phi_range),
                                     y=log10(v0_range),
                                     z=as.matrix(x[[2]]), 
                                     levels = 1, lwd=5, labels = "Cost of SM = 0.01",
                                     drawlabels = T, axes = FALSE, 
                                     frame.plot = FALSE, add = TRUE,labcex=1);
                 axis(1); axis(2); 
                 contour(x=log10(phi_range), 
                         y=log10(v0_range),
                         z=as.matrix(x[[1]]), 
                         levels = 1, lwd=2, labels = "Cost of SM = 0.001",
                         drawlabels = T, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,labcex=1);
                 contour(x=log10(phi_range),
                         y=log10(v0_range),
                         z=as.matrix(x[[3]]), 
                         levels = 1, lwd=2, labels = "Cost of SM = 0.1",
                         drawlabels = T, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,labcex=1);
                 #abline(v=1,lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 #abline(v=3,lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 #abline(h=log10(3e8),lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 #abline(h=log10(2e10),lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 points(1,log10(3e8),pch=21,bg="red",cex=2);
                 points(3,log10(3e10),pch=21,bg="red",cex=2);
                 text(x=-3,y=3,labels="CRISPR Favored",cex=1.5)
                 text(x=-3,y=11,labels="SM Favored",cex=1.5)},
               key.axes = "",labcex=1)
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()


