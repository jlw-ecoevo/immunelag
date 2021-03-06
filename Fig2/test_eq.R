

autoimmunitySystem <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    dR <- w*(r0-R) - e*v*R*(C+M)/(z+R)
    dC <- (v*R/(z+R) - w - delta*V)*C + phi*L
    dL <-  delta*C*V - phi*L - w*L
    dM <- ((1-k)*v*R/(z+R) - w)*M
    dV <- w*(v0-V)- delta*C*V
    
    # return the rate of change
    list(c(dR,dV,dC,dL,dM))
  }) # end with(as.list ...
}

testEq <- function(eq,parameters){
  times <- seq(0, 1e7, by = 1e6)
  state <- eq
  out <- ode(y = state, 
             times = times, 
             func = autoimmunitySystem, 
             parms = parameters,
             method="lsoda")
  return(tail(out,n=1)[,c("R","V","C","L","M")])
}