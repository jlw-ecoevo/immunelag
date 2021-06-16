
# JLW - 2020

# Infection dynamics under repeated outbreaks
# Panels D and E

library(deSolve)
library(ggsci)
library(scales)
library(wesanderson)
library(dplyr)


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
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D+Df)*V
    
    if(shape_return==2){
      dDf <- (v*R/(z+R) - w)*Df + mu*delta*V*U + phi*Id - zeta*Df*(Df/(Df+V))
      dD <-  (v*R/(z+R) - w - delta*V)*D + zeta*Df*(Df/(Df+V))
    } else if(shape_return==1){
      dDf <- (v*R/(z+R) - w)*Df + mu*delta*V*U + phi*Id - zeta*Df
      dD <-  (v*R/(z+R) - w - delta*V)*D + zeta*Df
    } else if(shape_return==0){
      dDf <- (v*R/(z+R) - w)*Df + mu*delta*V*U + phi*Id
      dD <-  (v*R/(z+R) - w - delta*V)*D
    } else {
      stop("Please choose valid shape for return from highly expressed state")
    }
    
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
                zeta=0,
                shape_return=2){
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
                  zeta=zeta,
                  shape_return=shape_return)
  return(parameters)
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
summaryIterateInducible <- function(time_to_outbreak,frac_suscept_outbreak,parameters,state){
  out <- iterateOutbreakInducible(time_to_outbreak = time_to_outbreak,
                              parameters = parameters,
                              state = state,
                              frac_suscept_outbreak = frac_suscept_outbreak,
                              n_iter=10)
  return(matrix(c(tail(out[,"D"],n=1)+tail(out[,"Id"],n=1)+tail(out[,"Df"],n=1),tail(out[,"SM"],n=1)),nrow=1))
}

## *** Takes a while to run

# parameter ranges to sweep over
t_range <- c(10^seq(0,4,.01))
f_range <- c(10^seq(-3,0,0.01))
sim_grid_A <- data.frame(Var1=t_range,Var2=0.5)
sim_grid_B <- data.frame(Var1=24,Var2=f_range)
sim_grid <- rbind(sim_grid_A,sim_grid_B)

# Solve w/ upregulation, nonlinear return
sim_grid_inducible <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e-1,shape_return=2)
D_mat_inducible <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible$D <- D_mat_inducible[1,]
sim_grid_inducible$SM <- D_mat_inducible[2,]

# Solve w/ upregulation, linear return
sim_grid_inducible_l <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e-1,shape_return=1)
D_mat_inducible <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible_l$D <- D_mat_inducible[1,]
sim_grid_inducible_l$SM <- D_mat_inducible[2,]

# Solve w/ upregulation, no return
sim_grid_inducible_nr <- sim_grid
state <- getInitialInducible(V=100)
parameters <- getParameters(v0=100,kappa=0.01,phi=10,zeta=1e-1,shape_return=0)
D_mat_inducible <- mapply(FUN=summaryIterateInducible,sim_grid_inducible$Var1,sim_grid_inducible$Var2, MoreArgs = list(parameters=parameters, state=state))
sim_grid_inducible_nr$D <- D_mat_inducible[1,]
sim_grid_inducible_nr$SM <- D_mat_inducible[2,]


# Save output
setwd("~/immunelag/Fig4")
save(sim_grid_inducible,
     sim_grid_inducible_l,
     sim_grid_inducible_nr,
     file="outbreak_iterated_compare_inducibility.RData")

setwd("~/immunelag/Fig4")
load("outbreak_iterated_compare_inducibility.RData")

 # Plot sweep over infected frac

sim_inducible <- sim_grid_inducible %>% subset(Var1==24)
sim_inducible_l <- sim_grid_inducible_l %>% subset(Var1==24)
sim_inducible_nr <- sim_grid_inducible_nr %>% subset(Var1==24)

Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible_l$D
Yf <- sim_inducible_l$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible_l <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible_nr$D
Yf <- sim_inducible_nr$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible_nr <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("compare_inducibility_rangefrac.pdf"),width=5,height=5)
plot(sim_inducible$Var2,rel_fit_inducible,
     type = "l",log="y",
     xlab = "Fraction of Host Susceptible to Outbreak",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#984ea3",
     ylim=c(1e-2,1e1),
     xlim=c(0.001,0.8))
lines(sim_inducible_nr$Var2,rel_fit_inducible_nr,lwd=3,col="#377eb8")
lines(sim_inducible_l$Var2,rel_fit_inducible_l,lwd=3,col="#e41a1c")
text(0.7,7e-1,"Linear Return",col="#e41a1c")
text(0.65,1.5e0,"Nonlinear Return",col="#984ea3")
text(0.62,3,"No Return",col="#377eb8")
text(0.25,8.5e-1,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()

# Plot sweep over time interval between outbreaks

sim_inducible <- sim_grid_inducible %>% subset(Var2==0.5)
sim_inducible_l <- sim_grid_inducible_l %>% subset(Var2==0.5)
sim_inducible_nr <- sim_grid_inducible_nr %>% subset(Var2==0.5)

Xf <- sim_inducible$D
Yf <- sim_inducible$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible_l$D
Yf <- sim_inducible_l$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible_l <- (Xprop*0.5)/((1-Xprop)*0.5)
Xf <- sim_inducible_nr$D
Yf <- sim_inducible_nr$SM
Xprop <- Xf/(Xf+Yf)
rel_fit_inducible_nr <- (Xprop*0.5)/((1-Xprop)*0.5)

setwd("~/immunelag/Fig4")
pdf(paste0("compare_inducibility_rangeinterval.pdf"),width=5,height=5)
plot(sim_inducible$Var1,rel_fit_inducible,
     type = "l",log="xy",
     ylim = c(1e-1,1e2),
     xlab = "Time Between Outbreaks (hours)",
     ylab = "Relative Fitness CRISPR x SM",lwd=3,
     col="#984ea3")
lines(sim_inducible_nr$Var1,rel_fit_inducible_nr,lwd=3,col="#377eb8")
lines(sim_inducible_l$Var1,rel_fit_inducible_l,lwd=3,col="#e41a1c")
text(5e2,9e0,"Linear Return",col="#e41a1c")
text(4.5e2,5,"Nonlinear Return",col="#984ea3")
text(2,2,"No Return",col="#377eb8")
text(160,1.2e0,"Rel. Fit. = 1",cex=0.6)
abline(h=1,lty=2)
dev.off()


