
# JLW - 2020
# Runs outbreak model w/ different growth rates

library(deSolve)
library(ggsci)
library(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(wesanderson)

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
    dD <- (v*R/(z+R) - w - delta*V)*D + mu*delta*U*V + phi*Id
    dId <-  delta*D*V - phi*Id - w*Id
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V

    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV))
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
    Id <- max(Id,0)
    dR <- w*(r0-R) - e*v*R*(U+D+SM)/(z+R)
    dU <- (v*R/(z+R) - delta*V - w)*U
    dI <- (1-mu)*delta*U*V - gam*I - w*I
    dD <- (v*R/(z+R) - w)*D + mu*delta*U*V
    dId <-  0
    dSM <- ((1-kappa)*v*R/(z+R) - w)*SM
    dV <- w*(v0-V) + beta*gam*I - delta*(U+D)*V

    # return the rate of change
    list(c(dR,dU,dI,dD,dId,dSM,dV))
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
                phi=100,
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


getCumRelFit <- function(v,phi){
  times <- seq(0, 336, by = .1)
  parameters <- getParameters(v=v,w=0.3,v0=100,kappa=0.01,phi=phi,r0=350)
  state <- getInitial(V=100)
  out <- ode(y = state,
             times = times,
             func = autoimmunitySystem,
             parms = parameters,
             method="lsoda")
  fC <- (out[,"D"]+out[,"Id"])/(out[,"D"]+out[,"Id"]+out[,"SM"])
  fCnow <- fC[-1]
  fCbefore <- fC[1]
  relfit <- (fCnow*(1-fCbefore))/(fCbefore*(1-fCnow))

  return(data.frame(Time=out[-1,"time"],
                    RelFit=relfit,
                    v=v))
}

getCumRelFitNoLag <- function(v,phi){
  times <- seq(0, 336, by = .1)
  parameters <- getParameters(v=v,w=0.3,v0=100,kappa=0.01,phi=phi,r0=350)
  state <- getInitial(V=100)
  out <- ode(y = state,
             times = times,
             func = autoimmunitySystemNoLag,
             parms = parameters,
             method="lsoda")
  fC <- (out[,"D"]+out[,"Id"])/(out[,"D"]+out[,"Id"]+out[,"SM"])
  fCnow <- fC[-1]
  fCbefore <- fC[1]
  relfit <- (fCnow*(1-fCbefore))/(fCbefore*(1-fCnow))

  return(data.frame(Time=out[-1,"time"],
                    RelFit=relfit,
                    v=v))
}

# Cumulative RelFit

v_range <- c(0.4,0.8,1.2,2)
v_df <- data.frame()
for(i in 1:length(v_range)){
 v_df <- rbind(v_df,getCumRelFit(v_range[i],phi=10))
}

v_range <- c(0.4,0.8,1.2,2)
vnl_df <- data.frame()
for(i in 1:length(v_range)){
  vnl_df <- rbind(vnl_df,getCumRelFitNoLag(v_range[i],phi=10))
}

p1 <- ggplot(v_df,aes(x=Time,y=RelFit,color=factor(v))) +
  geom_line(size=1.5) +
  scale_y_log10(limit=c(1e-10,1e3)) +
  theme_base() +
  geom_hline(yintercept=1,lty=2) +
  labs(color="Max. Growth Rate (v)") +
  ylab("Relative Fitness CRISPR x SM") +
  xlab("Time (Hours)") +
  scale_color_manual(values = wes_palette("Zissou1")) +
  theme(legend.position = "none",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("With Lag")

p2 <- ggplot(vnl_df,aes(x=Time,y=RelFit,color=factor(v))) +
  geom_line(size=1.5) +
  scale_y_log10(limit=c(1e0,1e1)) +
  theme_base() +
  geom_hline(yintercept=1,lty=2) +
  labs(color="Max. Growth Rate (v)") +
  ylab("Relative Fitness CRISPR x SM") +
  xlab("Time (Hours)") +
  scale_color_manual(values = wes_palette("Zissou1")) +
  theme(legend.position = c(0.35,0.7),
        text = element_text(size=14),
        legend.background = element_rect(size = 0.5, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Lag")

setwd("~/immunelag/Fig5")
pdf("lag_growth_rate_fitness_longlag.pdf",width=8,height=4)
ggarrange(p1,
          p2,
          ncol=2,
          labels=c("(a)","(b)"))
dev.off()





# Cumulative RelFit Short Lag

v_range <- c(0.4,0.8,1.2,2)
v_df <- data.frame()
for(i in 1:length(v_range)){
  v_df <- rbind(v_df,getCumRelFit(v_range[i],phi=1000))
}

v_range <- c(0.4,0.8,1.2,2)
vnl_df <- data.frame()
for(i in 1:length(v_range)){
  vnl_df <- rbind(vnl_df,getCumRelFitNoLag(v_range[i],phi=1000))
}

p1 <- ggplot(v_df,aes(x=Time,y=RelFit,color=factor(v))) +
  geom_line(size=1.5) +
  scale_y_log10(limit=c(1e-8,1e3)) +
  theme_base() +
  geom_hline(yintercept=1,lty=2) +
  labs(color="Max. Growth Rate (v)") +
  ylab("Relative Fitness CRISPR x SM") +
  xlab("Time (Hours)") +
  scale_color_manual(values = wes_palette("Zissou1")) +
  theme(legend.position = "none",
        text = element_text(size=14),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("With Lag") 

p2 <- ggplot(vnl_df,aes(x=Time,y=RelFit,color=factor(v))) +
  geom_line(size=1.5) +
  scale_y_log10(limit=c(1e0,1e1)) +
  theme_base() +
  geom_hline(yintercept=1,lty=2) +
  labs(color="Max. Growth Rate (v)") +
  ylab("Relative Fitness CRISPR x SM") +
  xlab("Time (Hours)") +
  scale_color_manual(values = wes_palette("Zissou1")) +
  theme(legend.position = c(0.35,0.7),
        text = element_text(size=14),
        legend.background = element_rect(size = 0.5, colour = 1),
        plot.title = element_text(hjust = 0.5)) +
  ggtitle("No Lag")

setwd("~/immunelag/Fig5")
pdf("lag_growth_rate_fitness_shortlag.pdf",width=8,height=4)
ggarrange(p1,
          p2,
          ncol=2,
          labels=c("(a)","(b)"))
dev.off()
