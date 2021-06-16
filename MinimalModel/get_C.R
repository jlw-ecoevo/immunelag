
Complex <- function(real,imaginary) {complex(real=real,imaginary=imaginary)}
Power <- function(a,b) {a^b}
Sqrt <- function(x) {sqrt(as.complex(x))}


getC <- function(v0,phi){
  C1 <- (-500000*
           (-2559 - 8530*phi + 450*v0 + 
              Sqrt(6456681 + 43044540*phi + 71740900*Power(phi,2) - 
                     2303100*v0 - 7677000*phi*v0 + 202500*Power(v0,2))))/
    (3 + 10*phi)
  C2 <- (500000*(2559 + 8530*phi - 450*v0 + 
                   Sqrt(6456681 + 43044540*phi + 71740900*Power(phi,2) - 
                          2303100*v0 - 7677000*phi*v0 + 202500*Power(v0,2))))/
    (3 + 10*phi)
  return(c(C1,C2))
}