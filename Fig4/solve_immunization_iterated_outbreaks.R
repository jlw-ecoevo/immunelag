
# JLW - 2020

# Infection dynamics under repeated outbreaks
# Panels D and E

library(deSolve)
library(ggsci)
library(scales)
library(wesanderson)
library(dplyr)

# ODEs
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

# ODEs w/out lag
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
    dU <- ((1-2*mu/alpha)*v*R/(z+R) - delta*V - w)*U
    dI <- (1-mu)*delta*U*V - gam*I - w*I 
    dD <- ((1-2*mu/alpha)*v*R/(z+R) - w)*D + mu*delta*V*U 
    dId <-  0
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV))
  }) # end with(as.list ...
}

# ODEs for system w/ CRISPR that can be upregulated
autoimmunitySystemInducible <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    Df <- max(Df,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM+Df)/(z+R)
    dU <- (v*R/(z+R) - delta*V - w)*U
    dI <- (1-mu)*delta*U*V - gam*I - w*I
    dD <-  (v*R/(z+R) - w - delta*V)*D + zeta*Df
    dDf <- (v*R/(z+R) - w)*Df + mu*delta*V*U + phi*Id - zeta*Df
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D+Df)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dD,dDf,dId,dSM,dV))
  }) # end with(as.list ...
}


# Helper function
getParameters <- function(w=0.3,
                r0=350,
                e=5e-7,
                v=2,
                z=1,
                delta=1e-7,
                mu=1e-7,
                alpha=0.4,
                gam=4/3,
                beta=80,
                v0=0,
                phi=4/3,
                xi=1/6,
                kappa=0.05,
                zeta=0){
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
                  kappa=kappa,
                  zeta=zeta)
  return(parameters)
}

# Helper function
getInitial <- function(R=350,
                       U=0,
                       I=0,
                       D=5e7,
                       Id=0,
                       SM=5e7,
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

# Helper function
getInitialInducible <- function(R=350,
                       U=0,
                       I=0,
                       D=0,
                       Df=5e7,
                       Id=0,
                       SM=5e7,
                       V=0){
  state <- c(R=R,
             U=U,
             I=I,
             D=D,
             Df=Df,
             Id=Id,
             SM=SM,
             V=V)
  return(state)
}

# Run repeated outbreaks (affecting CRISPR and SM) at set intervals w/ set fraction
# of defended host converted to susceptible
iterateOutbreak <- function(time_to_outbreak,parameters,state,frac_suscept_outbreak,n_iter=5,solver="lsoda"){
  times <- seq(0, time_to_outbreak, by = 1)
  out <- numeric()
  for(i in 1:n_iter){
    out_i <- ode(y = state, 
               times = times, 
               func = autoimmunitySystem, 
               parms = parameters,
               method=solver)
    out_i[,"time"] <- out_i[,"time"]+((i-1)*time_to_outbreak)
    out <- rbind(out,out_i)
    state <- getInitial(R=tail(out_i[,"R"],n=1),
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak+tail(out_i[,"U"],n=1),
                        I=tail(out_i[,"I"],n=1),
                        D=(1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=(1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),
                        V=tail(out_i[,"V"],n=1))
  }
  return(out)
}
iterateOutbreakNoLag <- function(time_to_outbreak,parameters,state,frac_suscept_outbreak,n_iter=5,solver="lsoda"){
  times <- seq(0, time_to_outbreak, by = 1)
  out <- numeric()
  for(i in 1:n_iter){
    out_i <- ode(y = state, 
                 times = times, 
                 func = autoimmunitySystemNoLag, 
                 parms = parameters,
                 method=solver)
    out_i[,"time"] <- out_i[,"time"]+((i-1)*time_to_outbreak)
    out <- rbind(out,out_i)
    state <- getInitial(R=tail(out_i[,"R"],n=1),
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak+tail(out_i[,"U"],n=1),
                        I=tail(out_i[,"I"],n=1),
                        D=(1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=(1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),
                        V=tail(out_i[,"V"],n=1))
  }
  return(out)
}
iterateOutbreakInducible <- function(time_to_outbreak,parameters,state,frac_suscept_outbreak,n_iter=5,solver="lsoda"){
  times <- seq(0, time_to_outbreak, by = 1)
  out <- numeric()
  for(i in 1:n_iter){
    out_i <- ode(y = state, 
                 times = times, 
                 func = autoimmunitySystemInducible, 
                 parms = parameters,
                 method=solver)
    out_i[,"time"] <- out_i[,"time"]+((i-1)*time_to_outbreak)
    out <- rbind(out,out_i)
    state <- getInitialInducible(R=tail(out_i[,"R"],n=1),
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"Df"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak+tail(out_i[,"U"],n=1),
                        I=tail(out_i[,"I"],n=1),
                        D=(1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),
                        Df=(1-frac_suscept_outbreak)*tail(out_i[,"Df"],n=1),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=(1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),
                        V=tail(out_i[,"V"],n=1))
  }
  return(out)
}

# Setup iterate outbreak and pull out final result
summaryIterate <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state){
  out <- iterateOutbreak(time_to_outbreak = time_to_outbreak,
                         parameters = parameters,
                         state = state,
                         frac_suscept_outbreak = frac_suscept_outbreak,
                         n_iter=10)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1),tail(out[,"SM"],n=1)),nrow=1))
}
summaryIterateNoLag <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state){
  out <- iterateOutbreakNoLag(time_to_outbreak = time_to_outbreak,
                         parameters = parameters,
                         state = state,
                         frac_suscept_outbreak = frac_suscept_outbreak,
                         n_iter=10)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1),tail(out[,"SM"],n=1)),nrow=1))
}
summaryIterateInducible <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state){
  out <- iterateOutbreakInducible(time_to_outbreak = time_to_outbreak,
                              parameters = parameters,
                              state = state,
                              frac_suscept_outbreak = frac_suscept_outbreak,
                              n_iter=10)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1)+tail(out[,"Df"],n=1),tail(out[,"SM"],n=1)),nrow=1))
}

## *** Takes a while to run

# initial conditions
state <- getInitial(V=100)

# parameter ranges to sweep over
t_range <- c(10^seq(0,4,.01))
f_range <- c(10^seq(-3,0,0.01))
sim_grid_A <- data.frame(Var1=t_range,Var2=0.5)
sim_grid_B <- data.frame(Var1=24,Var2=f_range)
sim_grid <- rbind(sim_grid_A,sim_grid_B)

# Solve w/ short lag
sim_grid_short <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=1000)
D_mat_short <- mapply(FUN=summaryIterate,sim_grid_short$Var1,sim_grid_short$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_short$D <- D_mat_short[1,]
sim_grid_short$SM <- D_mat_short[2,]

# Solve w/ long lag
sim_grid_long <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=10)
D_mat_long <- mapply(FUN=summaryIterate,sim_grid_long$Var1,sim_grid_long$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_long$D <- D_mat_long[1,]
sim_grid_long$SM <- D_mat_long[2,]

# Solve w/ no lag
sim_grid_none <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=0)
D_mat_none <- mapply(FUN=summaryIterateNoLag,sim_grid_none$Var1,sim_grid_none$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_none$D <- D_mat_none[1,]
sim_grid_none$SM <- D_mat_none[2,]

# Solve w/ upregulation
sim_grid_inducible <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e-1)
D_mat_inducible <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible$D <- D_mat_inducible[1,]
sim_grid_inducible$SM <- D_mat_inducible[2,]

# Solve w/ upregulation that decays quickly
sim_grid_inducible_short <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e1)
D_mat_inducible_short <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible_short$D <- D_mat_inducible_short[1,]
sim_grid_inducible_short$SM <- D_mat_inducible_short[2,]


# Save output
setwd("~/immunelag/Fig4")
save(sim_grid_short,
     sim_grid_long,
     sim_grid_none,
     sim_grid_inducible,
     sim_grid_inducible_short,
     file="outbreak_iterated.RData")

setwd("~/immunelag/Fig4")
load("outbreak_iterated.RData")

 # Plot sweep over infected frac

sim_long <- sim_grid_long %>% subset(Var1==24)
sim_short <- sim_grid_short %>% subset(Var1==24)
sim_none <- sim_grid_none %>% subset(Var1==24)
sim_inducible <- sim_grid_inducible %>% subset(Var1==24)

Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("outbreak_iterated_rangefrac.pdf"),width=5,height=5)
plot(sim_none$Var2,rel_fit_none,
     type = "l",log="y",
     xlab = "Fraction of Host Susceptible to Outbreak",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a",
     ylim=c(1e-2,1e1),
     xlim=c(0.001,0.8))
lines(sim_short$Var2,rel_fit_short,lwd=3,col="#377eb8")
lines(sim_long$Var2,rel_fit_long,lwd=3,col="#e41a1c")
lines(sim_inducible$Var2,rel_fit_inducible,lwd=3,col="#984ea3")
text(0.75,3e0,"No Lag",col="#4daf4a")
text(0.13,1e-1,"Long Lag",col="#e41a1c")
text(0.7,1.4e0,"Upregulation",col="#984ea3")
text(0.62,1e-1,"Short Lag",col="#377eb8")
text(0.25,8.5e-1,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()

# Plot sweep over time interval between outbreaks

sim_long <- sim_grid_long %>% subset(Var2==0.5)
sim_short <- sim_grid_short %>% subset(Var2==0.5)
sim_none <- sim_grid_none %>% subset(Var2==0.5)
sim_inducible <- sim_grid_inducible %>% subset(Var2==0.5)

Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("outbreak_iterated_rangeinterval.pdf"),width=5,height=5)
plot(sim_none$Var1,rel_fit_none,
     type = "l",log="xy",
     ylim = c(1e-1,1e2),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a")
lines(sim_none$Var1,rel_fit_short,lwd=3,col="#377eb8")
lines(sim_none$Var1,rel_fit_long,lwd=3,col="#e41a1c")
lines(sim_inducible$Var1,rel_fit_inducible,lwd=3,col="#984ea3")
text(1.3e1,3,"No Lag",col="#4daf4a")
text(5e3,9e0,"Long Lag",col="#e41a1c")
text(300,5,"Upregulation",col="#984ea3")
text(2,2.5e-1,"Short Lag",col="#377eb8")
text(160,1.2e0,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()

########### Plot output when upregulation decays quickly


sim_long <- sim_grid_long %>% subset(Var1==24)
sim_short <- sim_grid_short %>% subset(Var1==24)
sim_none <- sim_grid_none %>% subset(Var1==24)
sim_inducible <- sim_grid_inducible_short %>% subset(Var1==24)

Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("outbreak_iterated_rangefrac_shortinduction.pdf"),width=5,height=5)
plot(sim_none$Var2,rel_fit_none,
     type = "l",log="y",
     xlab = "Fraction of Host Susceptible to Outbreak",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a",
     ylim=c(1e-2,1e1),
     xlim=c(0.001,0.8))
lines(sim_short$Var2,rel_fit_short,lwd=3,col="#377eb8")
lines(sim_long$Var2,rel_fit_long,lwd=3,col="#e41a1c")
lines(sim_inducible$Var2,rel_fit_inducible,lwd=3,col="#984ea3")
text(0.75,3e0,"No Lag",col="#4daf4a")
text(0.13,1e-2,"Long Lag",col="#e41a1c")
text(0.28,2e-1,"Upregulation",col="#984ea3")
text(0.62,1e-1,"Short Lag",col="#377eb8")
text(0.25,8.5e-1,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()


sim_long <- sim_grid_long %>% subset(Var2==0.5)
sim_short <- sim_grid_short %>% subset(Var2==0.5)
sim_none <- sim_grid_none %>% subset(Var2==0.5)
sim_inducible <- sim_grid_inducible_short %>% subset(Var2==0.5)

Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("outbreak_iterated_rangeinterval_shortinduction.pdf"),width=5,height=5)
plot(sim_none$Var1,rel_fit_none,
     type = "l",log="xy",
     ylim = c(1e-1,1e2),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a")
lines(sim_none$Var1,rel_fit_short,lwd=3,col="#377eb8")
lines(sim_none$Var1,rel_fit_long,lwd=3,col="#e41a1c")
lines(sim_inducible$Var1,rel_fit_inducible,lwd=3,col="#984ea3")
text(1.3e1,3,"No Lag",col="#4daf4a")
text(5e3,9e0,"Long Lag",col="#e41a1c")
text(40,1e-1,"Upregulation",col="#984ea3")
text(2,2.5e-1,"Short Lag",col="#377eb8")
text(5e3,1.2e0,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()

