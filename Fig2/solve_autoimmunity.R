
# JLW - 2020

# Panel A for Figure 2

# Load packages
library(gplots)
library(dplyr)
library(ggplot2)
library(reshape2)
library(wesanderson)

# Parameter ranges to sweep over
alpha_range <- seq(0.1,1,0.01)
mu_range <- 10^seq(-10,0.5,0.05)
alpha_mu_combos <- expand.grid(alpha_range,mu_range)
names(alpha_mu_combos) <- c("alpha","mu")

# Invasion Condition
alpha_mu_combos$k <- 2*alpha_mu_combos$mu/alpha_mu_combos$alpha

# Reshape into matrix for plotting
k_mat <- reshape(alpha_mu_combos, idvar = "alpha", timevar = "mu", direction = "wide")

#Make a nice contour plot
setwd("~/immunelag/Fig2")
pdf(paste0("InvasionAutoimmunity_contour.pdf"),width=8,height=5)
par(mar=c(5.1, 5, 1, 0.2))
filled.contour(x=alpha_range,
               y=log10(mu_range),
               z=log10(as.matrix(k_mat[,-1])),
               xlab=expression("Maximum Rate of Genome Replication (" * alpha * ")"),
               ylab=expression("Spacer Acquisition Rate (" * log[10](mu) * ")"),
               levels=c(-10,-2,10),
               col = wes_palette("Moonrise3",2),
               labcex=1,
               cex.lab=1.5,
               cex.axis=1.5,
               plot.axes = { contour(x=alpha_range,
                                     y=log10(mu_range),
                                     z=log10(as.matrix(k_mat[,-1])),
                                     levels = -2, lwd=5, labels = "Cost of SM = 0.01",
                                     drawlabels = T, axes = FALSE, 
                                     frame.plot = FALSE, add = TRUE,labcex=1);
                 contour(x=alpha_range,
                         y=log10(mu_range),
                         z=log10(as.matrix(k_mat[,-1])),
                         levels = -3, lwd=2, labels = "Cost of SM = 0.001",
                         drawlabels = T, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,labcex=1);
                 contour(x=alpha_range,
                         y=log10(mu_range),
                         z=log10(as.matrix(k_mat[,-1])),
                         levels = -4, lwd=2, labels = "Cost of SM = 0.0001",
                         drawlabels = T, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,labcex=1);
                 contour(x=alpha_range,
                         y=log10(mu_range),
                         z=log10(as.matrix(k_mat[,-1])),
                         levels = -1, lwd=2, labels = "Cost of SM = 0.1",
                         drawlabels = T, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,labcex=1);
                 
                 abline(h=log10(4e-7),lty=2,col=wes_palette("Moonrise3",5)[3],lwd=1)
                 abline(h=log10(4e-6),lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 abline(h=log10(2e-5),lty=2,col=wes_palette("Moonrise3",5)[3],lwd=1)
                 abline(v=0.4,lty=2,col=wes_palette("Moonrise3",5)[3],lwd=4)
                 
                 points(0.4,log10(4e-6),pch=21,bg="red",cex=2);

                 text(x=0.7,y=-8,labels="CRISPR Favored",cex=1.5)
                 text(x=0.7,y=0,labels="SM Favored",cex=1.5);
                 axis(2); 
                 axis(1);},
               key.axes = "")
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()

