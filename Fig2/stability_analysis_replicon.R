library(deSolve)
library(wesanderson)

setwd("~/immunelag/Fig2")
source("get_R_replicon.R")
source("test_eq_replicon.R")

getJacobian <- function(state,parameters){
  with(as.list(c(state, parameters)),{
    #print(state)
    J <- matrix(0,nrow=5,ncol=5)
    J[1,1] <- -w - e*v*z*(M+C)/((z+R)^2)
    J[1,3] <- -e*v*R/(z+R)
    J[1,5] <- -e*v*R/(z+R)
    J[2,2] <- -w - delta*C
    J[2,3] <- -delta*V
    J[3,1] <- C*((alpha-2*mu)/alpha)*v*z/((z+R)^2)
    J[3,2] <- -delta*C
    J[3,3] <- ((alpha-2*mu)/alpha)*v*R/(z+R) - w - delta*V
    J[3,4] <- phi
    J[4,2] <- delta*C
    J[4,3] <- delta*V
    J[4,4] <- -phi - w
    J[5,1] <- M*(1-k)*v*z/((z+R)^2)
    J[5,5] <- (1-k)*v*z*R/(z+R) - w
    #print(J)
    return(J)
  })
}

getEqBoth <- function(parameters){
  with(as.list(c(parameters)),{
    Req <- z*w/(v-k*v-w)
    Veq <- (((2*mu-alpha)/alpha)*v*Req/(z+Req) + w)*((w+phi)/(-delta*w))
    Ceq <- w*(v0-Veq)/(Veq*delta)
    Leq <- delta*Veq*Ceq/(phi+w)
    Meq <- w*(r0-Req)*(z+Req)/(e*v*Req) - Ceq
    return(c(R = Req,
             V = Veq,
             C = Ceq,
             L = Leq,
             M = Meq))
  })
}

getEqM <- function(parameters){
  with(as.list(c(parameters)),{
    Req <- z*w/(v-k*v-w)
    Veq <- v0
    Ceq <- 0
    Leq <- 0
    Meq <- w*(r0-Req)*(z+Req)/(e*v*Req)
    return(c(R = Req,
             V = Veq,
             C = Ceq,
             L = Leq,
             M = Meq))
  })
}

getEqC <- function(parameters){
  with(as.list(c(parameters)),{
    Req <- getR(v0,phi)
    Ceq <- w*(r0-Req)*(z+Req)/(e*v*Req)
    Veq <- w*v0/(w+delta*Ceq)
    Leq <- delta*Veq*Ceq/(phi+w)
    Meq <- c(0,0,0)
    return(list(eq1 = c(R = Req[1],
                        V = Veq[1],
                        C = Ceq[1],
                        L = Leq[1],
                        M = Meq[1]),
                eq2 = c(R = Req[2],
                        V = Veq[2],
                        C = Ceq[2],
                        L = Leq[2],
                        M = Meq[2]),
                eq3 = c(R = Req[3],
                        V = Veq[3],
                        C = Ceq[3],
                        L = Leq[3],
                        M = Meq[3])))
  })
}

getEqNeither <- function(parameters){
  with(as.list(c(parameters)),{
    Req <- r0
    Veq <- v0
    Ceq <- 0
    Leq <- 0
    Meq <- 0
    return(c(R = Req,
             V = Veq,
             C = Ceq,
             L = Leq,
             M = Meq))
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
                  v0=v0,
                  phi=phi,
                  k=0.01)
  
  if(eq_type == "Both"){
    eq <- getEqBoth(parameters)
  } else if(eq_type == "M") {
    eq <- getEqM(parameters)
  } else if(eq_type == "C") {
    #warning("For this calculation all parameters assumed to be model defaultss except for v0 and phi, regardless of what you input into this function (formulas for roots pre-simplified in Mathematica)")
    eq <- getEqC(parameters)
  } else if(eq_type == "Neither") {
    eq <- getEqNeither(parameters)
  }else {
    stop("Please pick one of three EQ types (Both, M, C)")
  }
  
  if(eq_type != "C"){
    eq <- Re(eq)
    if(sum(abs(eq-testEq(eq,parameters)))>1e-3){
      warning("Computed equilibrium seems to be incorrect? (different from numerical system solved at large t starting at eq)")
      # print(sum(abs(eq-testEq(eq,parameters))))
    }
    J <- getJacobian(eq,parameters)
    return(getStability(eq,J))
    
  } else {
    if(sum(abs(eq[[1]]-testEq(eq[[1]],parameters)))>1e-3){
      warning("Computed equilibrium 1 seems to be incorrect? (different from numerical system solved at large t starting at eq)")
      #print(sum(abs(eq[[1]]-testEq(eq[[1]],parameters))))
    } else if(sum(abs(eq[[2]]-testEq(eq[[2]],parameters)))>1e-3){
      warning("Computed equilibrium 2 seems to be incorrect? (different from numerical system solved at large t starting at eq)")
      #print(sum(abs(eq[[2]]-testEq(eq[[2]],parameters))))
    } else if(sum(abs(eq[[3]]-testEq(eq[[3]],parameters)))>1e-3){
      warning("Computed equilibrium 3 seems to be incorrect? (different from numerical system solved at large t starting at eq)")
      #print(sum(abs(eq[[3]]-testEq(eq[[3]],parameters))))
    }
    eq <- lapply(eq,Re)
    J <- lapply(eq,getJacobian,parameters=parameters)
    return(mapply(getStability,eq,J))
  }

}


#Test
eqStability(1,1)
eqStability(1,1,eq_type="M")
eqStability(1,1,eq_type="C")
eqStability(1,1,eq_type="Neither")


#Multiple stable states in some cases
eqStability(1e6,0.1,eq_type="M")
eqStability(1e6,0.1,eq_type="C")


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
pdf(paste0("ModelEquilibria_contour_replicon.pdf"),width=8,height=5)
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
