

#Early infection dynamics

library(deSolve)
library(ggsci)
library(scales)
library(wesanderson)
library(dplyr)

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

autoimmunitySystemInducible <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    V <- max(V,0)
    U <- max(U,0)
    D <- max(D,0)
    Df <- max(Df,0)
    SM <- max(SM,0)
    I <- max(I,0)
    Iv <- max(Iv,0)
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM+Df)/(z+R)
    dU <- ((1-2*mu/alpha)*v*R/(z+R) - delta*V - w)*U
    dI <- delta*U*V - gam*I - w*I - mu*Iv
    dIv <- delta*U*V + gam*log(beta)*Iv - (gam*exp(xi*log(beta))/xi)*I - mu*(Iv^2)/(I+1e-16)
    dD <-  ((1-2*mu/alpha)*v*R/(z+R) - w - delta*V)*D + zeta*Df*(Df/(Df+delta*V))
    dDf <- ((1-2*mu/alpha)*v*R/(z+R) - w)*Df + mu*Iv + phi*Id - zeta*Df*(Df/(Df+delta*V))
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D+Df)*V
    
    # return the rate of change
    list(c(dR,dU,dI,dIv,dD,dDf,dId,dSM,dV))
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

getInitial <- function(R=350,
                       U=0,
                       I=0,
                       Iv=0,
                       D=5e7,
                       Id=0,
                       SM=5e7,
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

getInitialInducible <- function(R=350,
                       U=0,
                       I=0,
                       Iv=0,
                       D=0,
                       Df=5e7,
                       Id=0,
                       SM=5e7,
                       V=0){
  state <- c(R=R,
             U=U,
             I=I,
             Iv=Iv,
             D=D,
             Df=Df,
             Id=Id,
             SM=SM,
             V=V)
  return(state)
}


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
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak,
                        I=tail(out_i[,"I"],n=1),
                        Iv=tail(out_i[,"Iv"],n=1),
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
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak,
                        I=tail(out_i[,"I"],n=1),
                        Iv=tail(out_i[,"Iv"],n=1),
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
                        U=(tail(out_i[,"D"],n=1)+tail(out_i[,"Id"],n=1)+tail(out_i[,"Df"],n=1)+tail(out_i[,"SM"],n=1))*frac_suscept_outbreak,
                        I=tail(out_i[,"I"],n=1),
                        Iv=tail(out_i[,"Iv"],n=1),
                        D=(1-frac_suscept_outbreak)*tail(out_i[,"D"],n=1),
                        Df=(1-frac_suscept_outbreak)*tail(out_i[,"Df"],n=1),
                        Id=(1-frac_suscept_outbreak)*tail(out_i[,"Id"],n=1),
                        SM=(1-frac_suscept_outbreak)*tail(out_i[,"SM"],n=1),
                        V=tail(out_i[,"V"],n=1))
  }
  return(out)
}


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



state <- getInitial(V=100)
t_range <- c(10^seq(0,4,.01))
f_range <- c(10^seq(-3,0,0.01))

sim_grid_A <- data.frame(Var1=t_range,Var2=0.5)
sim_grid_B <- data.frame(Var1=24,Var2=f_range)
sim_grid <- rbind(sim_grid_A,sim_grid_B)

sim_grid_short <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=1000)
D_mat_short <- mapply(FUN=summaryIterate,sim_grid_short$Var1,sim_grid_short$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_short$D <- D_mat_short[1,]
sim_grid_short$SM <- D_mat_short[2,]

sim_grid_long <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=10)
D_mat_long <- mapply(FUN=summaryIterate,sim_grid_long$Var1,sim_grid_long$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_long$D <- D_mat_long[1,]
sim_grid_long$SM <- D_mat_long[2,]

sim_grid_none <- sim_grid
parameters <- getParameters(v0=100,kappa=0.01,phi=0)
D_mat_none <- mapply(FUN=summaryIterateNoLag,sim_grid_none$Var1,sim_grid_none$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_none$D <- D_mat_none[1,]
sim_grid_none$SM <- D_mat_none[2,]

sim_grid_inducible <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e-1)
D_mat_inducible <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible$D <- D_mat_inducible[1,]
sim_grid_inducible$SM <- D_mat_inducible[2,]


sim_grid_inducible_short <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e1)
D_mat_inducible_short <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible_short$D <- D_mat_inducible_short[1,]
sim_grid_inducible_short$SM <- D_mat_inducible_short[2,]



setwd("~/immunelag/Fig4")
save(sim_grid_short,
     sim_grid_long,
     sim_grid_none,
     sim_grid_inducible,
     sim_grid_inducible_short,
     file="outbreak_iterated.RData")

setwd("~/immunelag/Fig4")
load("outbreak_iterated.RData")


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

########### Short Induction


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
text(0.28,5e-2,"Upregulation",col="#984ea3")
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

