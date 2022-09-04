###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 3.                                ##
## Spatial components and transport          ##
## Chapter 3.6.1.                             ##
## Autocatalysis in a flow-through           ##
## stirred tank                              ##
###############################################


# the model
autocatalysis <- function(t,state,pars)
{

with (as.list(c(state,pars)),

{
 dA <- dr*(Ain-A)-k*A*B
 dB <- dr*(Bin-B)+k*A*B
 dC <- -dr*C +k*A*B
 return (list(c(dA,dB,dC)))
})

}

# model application
times <- seq(0,300,1)
state <- c(A=1,B=1,C=0)
parms <- c(Ain =1,Bin = 0.1, k = 0.05, dr = 0.05)

require(deSolve)
out   <- as.data.frame(ode(state,times,autocatalysis,parms))
ylim  <- range(c(out$A,out$B,out$C))

plot(out$time,out$A,xlab="time",ylab="concentration",
      lwd=2,type="l",ylim=ylim,main="autocatalysis")
lines(out$time,out$B,lwd=2,lty=2)
lines(out$time,out$C,lwd=2,lty=3)

legend("topright",c("A","B","C"),lwd=2,lty=c(1,2,3))
