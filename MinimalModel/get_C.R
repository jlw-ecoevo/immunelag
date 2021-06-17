
Complex <- function(real,imaginary) {complex(real=real,imaginary=imaginary)}
Power <- function(a,b) {a^b}
Sqrt <- function(x) {sqrt(as.complex(x))}


getC <- function(v0,phi){
  C1 <- (5000*(254100 + 847000*phi - 
                 Sqrt(2)*Sqrt((3 + 10*phi)*
                                (10914135000 + 36380450000*phi - 9*v0))))/
    (3 + 10*phi)
  C2 <- (5000*(254100 + 847000*phi + 
                Sqrt(2)*Sqrt((3 + 10*phi)*
                               (10914135000 + 36380450000*phi - 9*v0))))/
  (3 + 10*phi)
  return(c(C1,C2))
}