###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.9.2. Project                    ##
## Fisheries model with Allee effect         ##
###############################################

require(shape)
require(rootSolve)

# parameter values
r   <- 0.05
K   <- 10    
K0  <- 1.

# the model : density-dependent growth, allee effect and linear mortality
rate <- function(x,mrt=0.02) r*x*(1-x/K)*(x/K0-1)-mrt*x

# plot the function for different values of mortality
xseq   <- seq(0,10,length=500)
mrtseq <- seq(0.0,0.1,by=0.02)
mat  <- outer (X=xseq,Y=mrtseq, function(X,Y) rate(x=X,mrt=Y))
matplot(xseq,mat,xlab="x",ylab="dx/dt",type="l",lty=1:10,col=1,main="fisheries model with Allee effect") 
abline(h=0,lty=2)                # add 0-axis
legend("bottomleft",legend=mrtseq,title="mrt",col=1,lty=1:10)
writelabel("C")

# bifurcation diagram
# open a plot
mrtseq <- seq(0.0,0.15,by=0.0005)
plot(0,xlim=range(mrtseq),ylim=c(0,10),type="n",
     xlab="fishing mortality",ylab="Equilibrium biomass", 
     main="fisheries model with Allee effect")

# loop over mortality values, find all roots in interval and 
# use sign of Jacobian to estimate if stable (neg), saddle(0) or unstable(pos)
for (mrt in mrtseq) {
  equi   <- uniroot.all(f=rate,interval=c(0,10),mrt=mrt)
  jac    <- diag(gradient(f =rate,x=equi,mrt=mrt))
  eig    <- sign(jac) 

  points(rep(mrt,length(equi)),equi,pch=22,
   col=c("darkgrey","black","lightgrey")[eig+2],
   bg =c("darkgrey","black","lightgrey")[eig+2]) 
}

legend("topright",pch=22,pt.cex=2,c("stable","unstable"),
col=c("darkgrey","lightgrey"),pt.bg=c("darkgrey","lightgrey"))
writelabel("D")