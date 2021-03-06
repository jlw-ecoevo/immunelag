
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


# Helper function
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



# Run repeated outbreaks (affecting CRISPR and SM) at set intervals w/ set fraction
# of defended host converted to susceptible
iterateOutbreak <- function(time_to_outbreak,parameters,state,frac_suscept_outbreak,n_iter=5,solver="lsoda"){
  times <- seq(0, time_to_outbreak, by = 1)
  out <- numeric()
  for(i in 1:n_iter){
    # print(state)
    out_i <- ode(y = state, 
               times = times, 
               func = autoimmunitySystem, 
               parms = parameters,
               method=solver)
    out_i[,"time"] <- out_i[,"time"]+((i-1)*time_to_outbreak)
    out <- rbind(out,out_i)
    state <- getInitial(R=tail(out_i[,"R"],n=1),
                        U=(tail(out_i[,"D"],n=1) + tail(out_i[,"Id"],n=1) + tail(out_i[,"SM"],n=1))*frac_suscept_outbreak + tail(out_i[,"U"],n=1),
                        I=tail(out_i[,"I"],n=1),
                        D=max((1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),0),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=max((1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),1),
                        V=as.numeric(parameters["v0"]))
  }
  return(out)
}
iterateOutbreakNoLag <- function(time_to_outbreak,parameters,state,frac_suscept_outbreak,n_iter=5,solver="lsoda"){
  times <- seq(0, time_to_outbreak, by = 1)
  out <- numeric()
  for(i in 1:n_iter){
    # print(state)
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
                        D=max((1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),0),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=max((1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),1),
                        V=as.numeric(parameters["v0"]))
  }
  return(out)
}


# Setup iterate outbreak and pull out final result
summaryIterate <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state,n_iter=10){
  out <- iterateOutbreak(time_to_outbreak = time_to_outbreak,
                         parameters = parameters,
                         state = state,
                         frac_suscept_outbreak = frac_suscept_outbreak,
                         n_iter=n_iter)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1),tail(out[,"SM"],n=1)),nrow=1))
  # out <- out[out[,"time"]>((n_iter-5)*time_to_outbreak),]
  # return(matrix(c(mean(out[,"D"]+out[,"Id"]),mean(out[,"SM"])),nrow=1))
}
summaryIterateNoLag <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state,n_iter=10){
  out <- iterateOutbreakNoLag(time_to_outbreak = time_to_outbreak,
                         parameters = parameters,
                         state = state,
                         frac_suscept_outbreak = frac_suscept_outbreak,
                         n_iter=n_iter)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1),tail(out[,"SM"],n=1)),nrow=1))
  # out <- out[out[,"time"]>((n_iter-5)*time_to_outbreak),]
  # return(matrix(c(mean(out[,"D"]+out[,"Id"]),mean(out[,"SM"])),nrow=1))
}

## *** Takes a while to run


# parameter ranges to sweep over
t_range <- c(10^seq(0,4,.05))
sim_grid_A <- data.frame(Var1=t_range,Var2=0.5)
sim_grid_B <- data.frame(Var1=t_range,Var2=0.99)
sim_grid_C <- data.frame(Var1=t_range,Var2=1)
sim_grid <- rbind(sim_grid_A,sim_grid_B,sim_grid_C)


# No Priming -----------------------------------

parameters <- getParameters(v0=0,kappa=0.01)
state <- getInitial(V=0,D=0,SM=1)
times <- seq(0, 24, by = .1)
out <- ode(y = state, 
    times = times, 
    func = autoimmunitySystemNoLag, 
    parms = parameters,
    method="lsoda")
plot(out[,"time"],out[,"SM"],type="l",lwd=2,log="y",ylim=c(1,1e9))
abline(v=12,lty=2)

# initial conditions
state <- getInitial(V=100)

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



# Save output
setwd("~/immunelag/Fig4")
save(sim_grid_short,
     sim_grid_long,
     sim_grid_none,
     file="outbreak_iterated_defective.RData")

setwd("~/immunelag/Fig4")
load("outbreak_iterated_defective.RData")


# Plot sweep over time interval between outbreaks ------------------------------

# Not primed, f=0.5
sim_long <- sim_grid_long %>% subset(Var2==0.5)
sim_short <- sim_grid_short %>% subset(Var2==0.5)
sim_none <- sim_grid_none %>% subset(Var2==0.5)
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

# Not primed, f=0.99
sim_long <- sim_grid_long %>% subset(Var2==0.99)
sim_short <- sim_grid_short %>% subset(Var2==0.99)
sim_none <- sim_grid_none %>% subset(Var2==0.99)
Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long_f99 <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short_f99 <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none_f99 <- (Xprop*0.5)/((1-Xprop)*0.5)

# Not primed, f=1
sim_long <- sim_grid_long %>% subset(Var2==1)
sim_short <- sim_grid_short %>% subset(Var2==1)
sim_none <- sim_grid_none %>% subset(Var2==1)
Xf <- sim_long$D
Yf <- sim_long$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_long_f1 <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_short$D
Yf <- sim_short$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_short_f1 <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_none$D
Yf <- sim_none$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_none_f1 <- (Xprop*0.5)/((1-Xprop)*0.5)



setwd("~/immunelag/Fig4")
pdf(paste0("outbreak_iterated_rangeinterval_defective.pdf"),width=15,height=5)
par(mfrow=c(1,3))

plot(sim_none$Var1,rel_fit_none,
     type = "l",log="xy",
     ylim = c(1e-1,1e5),
     xlim = c(15,1e4),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a",
     main="f=0.5")
lines(sim_none$Var1,rel_fit_short,lwd=3,col="#377eb8")
lines(sim_none$Var1,rel_fit_long,lwd=3,col="#e41a1c")
text(25,5,"No Lag",col="#4daf4a")
text(3.5e3,9e0,"Long Lag",col="#e41a1c")
text(90,1e3,"Short Lag",col="#377eb8")
text(160,1.2e0,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)

plot(sim_none$Var1,rel_fit_none_f99,
     type = "l",log="xy",
     ylim = c(1e-1,1e5),
     xlim = c(15,1e4),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a",
     main="f=0.99")
lines(sim_none$Var1,rel_fit_short_f99,lwd=3,col="#377eb8")
lines(sim_none$Var1,rel_fit_long_f99,lwd=3,col="#e41a1c")
abline(h=1,lty=2)

plot(sim_none$Var1,rel_fit_none_f1,
     type = "l",log="xy",
     ylim = c(1e-1,1e5),
     xlim = c(15,1e4),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#4daf4a",
     main="f=1")
lines(sim_none$Var1,rel_fit_short_f1,lwd=3,col="#377eb8")
lines(sim_none$Var1,rel_fit_long_f1,lwd=3,col="#e41a1c")
abline(h=1,lty=2)
abline(v=12,lty=2)

dev.off()

