


library(deSolve)
library(dplyr)
library(ggplot2)
library(ggpubr)


bootMeanCI <- function(x, n = 1000){
  boot.samples <- matrix(sample(x, size = n*length(x), replace = TRUE), 
                         n, 
                         length(x))
  boot.means <- apply(boot.samples, 1, mean)
  return(quantile(x,c(0.025,0.975)))
}


lagSystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dR <- - e*v*R*(X+Y)/(z+R)
    dX <- ((1-2*mu/alpha)*v*R/(z+R) - delta*V)*X + phi*Ix
    dY <- ((1-k)*v*R/(z+R))*Y
    dIx <- delta*X*V - phi*Ix
    dV <- -delta*X*V
    # return the rate of change
    list(c(dR,dX,dY,dIx,dV))
  }) # end with(as.list ...
}

competeXY <- function(V0,phi,dpi=1,re_add_phage=T,daylen=24,kappa=0){
  parameters <- c(e=5e-7,
                  v=2,
                  z=1,
                  delta=1e-7,
                  mu=4e-6,
                  alpha=0.4,
                  k=kappa,
                  delta=1e-7,
                  phi=phi)
  state <- c(R=500,
             X=max(1e7-V0,0),
             Y=1e7,
             Ix=min(V0,1e7),
             V=V0)
  times <- seq(0, daylen, by = 1)
  
  for(i in 1:dpi){
    out <- ode(y = state, 
               times = times, 
               func = lagSystem, 
               parms = parameters,
               method="lsoda")
    if(re_add_phage==T){
      state <- c(R=500+tail(out[,"R"],n=1)/100,
                 X=tail(out[,"X"],n=1)/100,
                 Y=tail(out[,"Y"],n=1)/100,
                 Ix=tail(out[,"Ix"],n=1)/100,
                 V=V0+tail(out[,"V"],n=1)/100)
    } else{
      state <- c(R=500+tail(out[,"R"],n=1)/100,
                 X=tail(out[,"X"],n=1)/100,
                 Y=tail(out[,"Y"],n=1)/100,
                 Ix=tail(out[,"Ix"],n=1)/100,
                 V=tail(out[,"V"],n=1)/100)
    }
  }
  
  Xf <- tail(out[,"X"],n=1)
  Yf <- tail(out[,"Y"],n=1)
  Xprop <- Xf/(Xf+Yf)
  rel_fit <- (Xprop*0.5)/((1-Xprop)*0.5)
  return(rel_fit)
}

titreProfile <- function(phi,kappa,v_range=10^seq(0,11,1)){
  fitness <- numeric()
  for(v in v_range){
    fitness <- c(fitness,competeXY(V0=v,
                                   phi=phi,
                                   dpi=1,
                                   re_add_phage=F,
                                   daylen=24,
                                   kappa=kappa))
  }
  return(fitness)
}

setwd("~/immunelag/Fig3")
cur_exp <- read.csv("experiment.csv") %>% subset(select=c(Titre,Fitness))
cur_exp$Source <- "Current"
EW_exp <- data.frame(Titre =c(1e7/6,1e8/6,3e8/6), 
                     Fitness = c(1.95,0.28,0.014),
                     Source="Westra et al.")
EA_exp <- read.csv("alseth.csv",stringsAsFactors = F)  %>% 
  subset(Treatment %in% c("Monoculture (BIM2 + SM only)","Monoculture (BIM2 + SM)")) %>%
  subset(select=c(Phage.titre,CRISPR..BIM2..fitness))
names(EA_exp) <- c("Titre","Fitness")
EA_exp$Source <- "Alseth et al."
EA_exp$Titre <- gsub("[\\^]","e",EA_exp$Titre) %>% 
  gsub(pattern="10e",replace="1e") %>%
  as.numeric()
EA_exp$Titre <- EA_exp$Titre/6

anova.test <- aov(Fitness~factor(Titre),data=cur_exp)
summary(anova.test)
TukeyHSD(anova.test)


v_range <- 10^seq(3,11,0.1)
x <- titreProfile(phi=1e3,kappa=0,v_range=v_range)
model_df <- data.frame(v=v_range,
                       fitness=x,
                       stringsAsFactors = F)
v_range <- 10^seq(3,11,0.1)
x <- titreProfile(phi=4/3,kappa=0,v_range=v_range)
model_df_low <- data.frame(v=v_range,
                           fitness=x,
                           stringsAsFactors = F)



mean_exp <- cur_exp %>% group_by(Titre) %>% 
  summarise(MeanFitness=mean(Fitness),
            LowerCI=bootMeanCI(Fitness,n=1e5)[1],
            UpperCI=bootMeanCI(Fitness,n=1e5)[2])
p1 <- ggplot() + 
  geom_line(data=model_df,aes(x=v_range,y=fitness),color="#d95f02",lwd=1) +
  #geom_line(data=model_df_low,aes(x=v_range,y=fitness),color="red",lwd=0.5,alpha=0.5) +
  geom_point(data = mean_exp,
             aes(x=Titre,y=MeanFitness))  +
  geom_errorbar(data = mean_exp,
                aes(x=Titre,y=MeanFitness,ymin = LowerCI, ymax = UpperCI),width=0.1) + 
  scale_x_log10(limits=c(5e3,1e11)) + 
  ylim(0.1,1.4) +
  theme_pubclean() + 
  geom_jitter(data=cur_exp,aes(x=Titre,y=Fitness),alpha=0.25,width=0.1) +
  xlab("Phage Titre (PFU/mL)") + 
  ylab("Relative Fitness CRISPR x SM") + 
  geom_line(data=data.frame(x=c(1e8,1e10),y=c(1.1,1.1)),aes(x=x,y=y),alpha=0.5) +  
  geom_line(data=data.frame(x=c(1e6,1e10),y=c(1.2,1.2)),aes(x=x,y=y),alpha=0.5) + 
  geom_line(data=data.frame(x=c(1e4,1e10),y=c(1.3,1.3)),aes(x=x,y=y),alpha=0.5)  +
  geom_text(data=data.frame(x=c(1e7),y=c(1.32)),aes(x=x,y=y),alpha=0.5,label="**",size=5) +
  geom_text(data=data.frame(x=c(1e8),y=c(1.22)),aes(x=x,y=y),alpha=0.5,label="***",size=5) + 
  geom_text(data=data.frame(x=c(1e9),y=c(1.12)),aes(x=x,y=y),alpha=0.5,label="****",size=5) + 
  scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1)) + 
  geom_text(data=data.frame(x=c(3e10),y=c(0.75)),aes(x=x,y=y),color="#d95f02",label=expression(phi * "=" *1000),size=7)

## Alseth et al data

v_range_EA <- 10^seq(0,10,0.1)
x <- titreProfile(phi=10,kappa=0,v_range=v_range_EA)
model_df_EA <- data.frame(v=v_range_EA,
                       fitness=x,
                       stringsAsFactors = F)

mean_exp_EA <- EA_exp %>% group_by(Titre) %>% 
  summarise(MeanFitness=mean(Fitness),
            LowerCI=bootMeanCI(Fitness,n=1e5)[1],
            UpperCI=bootMeanCI(Fitness,n=1e5)[2])

p2 <- ggplot() + 
  geom_line(data=model_df_EA,aes(x=v,y=fitness),color="#7570b3",lwd=1) +
  geom_point(data = mean_exp_EA,
             aes(x=Titre+1,y=MeanFitness))  +
  geom_errorbar(data = mean_exp_EA,
                aes(x=Titre+1,y=MeanFitness,ymin = LowerCI, ymax = UpperCI),width=0.1) + 
  scale_x_log10() +
  theme_pubclean() + 
  geom_jitter(data=EA_exp,aes(x=Titre+1,y=Fitness),alpha=0.25,width=0.1) +
  xlab("Phage Titre (PFU/mL)") + 
  ylab("Relative Fitness CRISPR x SM")  + 
  scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5),limit=c(0,2.75)) +
  ggtitle("Alseth et al. 2019 (Fig 2a)") + 
  geom_text(data=data.frame(x=c(1e8),y=c(1)),aes(x=x,y=y),color="#7570b3",label=expression(phi * "=" * 10),size=7)



## Westra et al data

v_range_EW <- 10^seq(0,10,0.1)
x <- titreProfile(phi=10,kappa=0.15,v_range=v_range_EW)
model_df_cost_EW <- data.frame(v=v_range_EW,
                           fitness=x,
                           stringsAsFactors = F)
v_range_EW <- 10^seq(0,10,0.1)
x <- titreProfile(phi=10,kappa=0,v_range=v_range_EW)
model_df_EW <- data.frame(v=v_range_EW,
                           fitness=x,
                           stringsAsFactors = F)


p3 <- ggplot() + 
  geom_line(data=model_df_cost_EW,aes(x=v,y=fitness),color="#7570b3",lwd=0.5,lty=2) +
  geom_line(data=model_df_EW,aes(x=v,y=fitness),color="#7570b3",lwd=1) +
  geom_point(data = data.frame(x = c(1e7,1e8,3e8), 
                               y = c(1.95,0.28,0.014),
                               CI = c(0.62,0.03,0.011)),
             aes(x=x,y=y))  +
  geom_errorbar(data = data.frame(x = c(1e7,1e8,3e8), 
                                  y = c(1.95,0.28,0.014),
                                  CI = c(0.62,0.03,0.011)),
                aes(x=x,y=y,ymin = y - CI, ymax = y + CI),width=0.1) + 
  scale_x_log10(limits=c(1e0,1e10)) + 
  theme_pubclean() + 
  xlab("Phage Titre (PFU/mL)") + 
  ylab("Relative Fitness CRISPR x SM")   + 
  scale_y_continuous(breaks=c(0,0.5,1,1.5,2,2.5),limit=c(0,2.75))+
  ggtitle("Westra et al. 2015 (Fig S3F)") + 
  geom_text(data=data.frame(x=c(1e1),y=c(1.9)),
          aes(x=x,y=y),alpha=1,label=expression(kappa * "=" * 0.15),
          size=3,color="black") +
  geom_text(data=data.frame(x=c(1e1),y=c(1.05)),
            aes(x=x,y=y),alpha=1,label=expression(kappa * "=" * 0),
            size=3,color="black") + 
  geom_text(data=data.frame(x=c(2e5),y=c(0.7)),aes(x=x,y=y),color="#7570b3",label=expression(phi * "=" * 10),size=7)

## Put it together


setwd("~/immunelag/Fig3")
pdf(paste0("CRISPRxSM_allexperiments.pdf"),width=6,height=10)
ggarrange(p1,
          ggarrange(p2,p3,ncol=2,labels=c("(b)","(c)")),
          nrow=2,
          labels=c("(a)",""))
dev.off()
