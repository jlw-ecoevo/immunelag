library(deSolve)
library(wesanderson)

setwd("~/immunelag/MinimalModel")

getJacobian <- function(state,parameters){
  with(as.list(c(state, parameters)),{
    #print(state)
    J <- matrix(0,nrow=4,ncol=4)
    J[1,1] <- r*(1 - ((2*C+M)/k)) - delta*V
    J[1,2] <- -r*C/k
    J[1,3] <- phi
    J[1,4] <- -delta*C
    J[2,1] <- -r*(1-kappa)*M/k
    J[2,2] <- r*(1-kappa)*(1-((C+2*M)/k))
    J[3,1] <- delta*V
    J[3,3] <- -phi
    J[3,4] <- delta*C
    J[4,1] <- -delta*V
    J[4,4] <- -delta*C-w
    #print(J)
    return(J)
  })
}


getEqC <- function(parameters){
  with(as.list(c(parameters)),{
    Veq <- w*v0/(delta*k - w)
    Ceq <- k
    Leq <- k*Veq/phi
    Meq <- 0
    return(c(C = Ceq,
             M = Meq,
             L = Leq,
             V = Veq))
  })
}

getEqM <- function(parameters){
  with(as.list(c(parameters)),{
    Veq <- v0
    Ceq <- 0
    Leq <- 0
    Meq <- k
    return(c(C = Ceq,
             M = Meq,
             L = Leq,
             V = Veq))
  })
}


reig <- function(x) {Re(eigen(x)$values)}
getStability <- function(eq,J){
  if(sum(reig(J)>=0)==0 & sum(eq<0)==0){#Stable eq
    return(1)
  } else if(sum(reig(J)>0)>0 & sum(eq<0)==0){#Unstable eq
    return(-1)
  } else if(sum(reig(J)==0)>0 & sum(reig(J)>0)==0 & sum(eq<0)==0){#Unclear
    return(2)
  }  else if(sum(eq<0)>0){#Not a real eq
    return(3)
  }
}

eqStability <- function(v0,phi,eq_type = "Both") {
  
  parameters <- c(r=1,
                  delta=1e-7,
                  v0=v0,
                  phi=phi,
                  kappa=0.01,
                  k=1e9,
                  w=0.3)
  
  if(eq_type == "M") {
    eq <- getEqM(parameters)
  } else if(eq_type == "C") {
    eq <- getEqC(parameters)
  } else {
    stop("Please pick one of two EQ types (M, C)")
  }
  # print(eq)
 
  J <- getJacobian(eq,parameters)
  print(J)
  print(reig(J))
  return(getStability(eq,J))
}


#Test
eqStability(1,1)
eqStability(0,1,eq_type="M")
eqStability(0,1,eq_type="C")


#Multiple stable states in some cases
eqStability(1e10,0.1,eq_type="M")
eqStability(1e10,0.1,eq_type="C")


# Phase Diagrams --------------------------------------------------------------



phi_range <- (10^seq(-2,5,0.1))
v0_range <- 10^seq(0,12,0.1)

pv_combos <- expand.grid(phi_range,v0_range)
names(pv_combos) <- c("phi","v0")
x <- mapply(eqStability,phi=pv_combos$phi,v0=pv_combos$v0)
table(x)
pv_combos$B <- x
B_mat <- reshape(pv_combos, idvar = "phi", timevar = "v0", direction = "wide")
table(unlist(B_mat[,-1])) # Never stable



phi_range <- (10^seq(-2,5,0.1))
v0_range <- 10^seq(0,12,0.1)

pv_combos <- expand.grid(phi_range,v0_range)
names(pv_combos) <- c("phi","v0")
x <- mapply(eqStability,phi=pv_combos$phi,v0=pv_combos$v0,MoreArgs = list(eq_type="M"))
table(x)
pv_combos$M <- x
M_mat <- reshape(pv_combos, idvar = "phi", timevar = "v0", direction = "wide")

pv_combos <- expand.grid(phi_range,v0_range)
names(pv_combos) <- c("phi","v0")
x <- mapply(eqStability,phi=pv_combos$phi,v0=pv_combos$v0,MoreArgs = list(eq_type="C"))
table(x[1,])
table(x[2,])
table(x[3,])
pv_combos$C <- x[3,]
C_mat <- reshape(pv_combos, idvar = "phi", timevar = "v0", direction = "wide")


eq_mat <- C_mat[,-1] == 1
eq_mat[C_mat[,-1] == 1 & M_mat[,-1] == 1] <- 2
eq_mat[C_mat[,-1] != 1 & M_mat[,-1] == 1] <- 3

eq_mat_rev <- eq_mat[nrow(eq_mat):1,]


x <- rev(1/phi_range)
xtick <- 10^(seq(-5,2,2))
x <- log(x)
y <- v0_range
ytick <- 10^(seq(0,12,3))
y <- log(y)
pdf(paste0("ModelEquilibria_contour_minimal.pdf"),width=8,height=5)
par(mar=c(5.1, 5, 1, 0.2))
filled.contour(x=x,
               y=y,
               z=eq_mat_rev,
               xlab=expression("Lag Length (" * 1/phi * ")"),
               ylab=expression("Envrionmental Viral Pool (" * v[0] * ")"),
               levels=c(0,1.5,2.5,4),
               col = wes_palette("Moonrise3",3),
               cex.lab=1.5,
               cex.axis=1.5,
               plot.axes = {axis(1, at=log(xtick), label=xtick);
                            axis(2, at=log(ytick), label=ytick)
                 contour(x=x,
                         y=y,
                         z=eq_mat_rev,
                         levels=c(0,1.5,2.5,4), 
                         lwd=1, 
                         drawlabels = F, axes = FALSE, 
                         frame.plot = FALSE, add = TRUE,
                         col="black")
                 text(x=2,y=7,labels="CRISPR Only",cex=1.25)
                 text(x=2,y=19,labels="SM Only",cex=1.25)
                 text(x=2,y=13,labels="CRISPR or SM",cex=1.25)
                 text(x=-3,y=23.6,labels=expression(v[0] * "=" * 10^9),cex=0.8)
                 abline(h=log(1e9),lty=2,lwd=2)})
par(mar=c(5.1, 4.1, 4.1, 2.1))
dev.off()
