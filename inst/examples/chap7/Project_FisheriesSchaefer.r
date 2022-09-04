###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.9.1. Project                    ##
## The Schaefer model of fisheries           ##
###############################################

require(shape)
windows()
par(mfrow=c(2,2))
# parameter values
r   <- 0.05
K   <- 10    
mrt <- 0.02   

# the model : density-dependent growth and sigmoid-type mortality rate
rate <- function(x,mrt=0.02) r*x*(1-x/K)-mrt*x

# plot the function for different values of mortality
xseq   <- seq(0,10,length=500)
mrtseq <- seq(0.0,0.1,by=0.02)
mat  <- outer (X=xseq,Y=mrtseq, function(X,Y) rate(x=X,mrt=Y))
matplot(xseq,mat,xlab="x",ylab="dx/dt",type="l",lty=1:10,col=1,main="Schaefer model") 
abline(h=0,lty=2)                # add 0-axis
legend("bottomleft",legend=mrtseq,title="mrt",col=1,lty=1:10)
writelabel("A")

#windows()
par(mar=c(5.1,4.1,4.1,4.1))
mrtseq <- seq(0.0,0.05,by=0.0005)
equi  <- K*(1-mrtseq/r)
Yield <- equi*mrtseq
plot(mrtseq,equi,type="l",lwd=2,xlab="fishing mortality",ylab="Equilibrium biomass", 
     main="Schaefer model")
par(new=TRUE)
plot(mrtseq,Yield,axes=FALSE,xlab="",ylab="",type="l",lwd=1)
axis(4)
mtext(side=4,line=2.5,"Fisheries yield",cex=0.8)
mmax <-mrtseq[which.max(Yield)]
xmax <- equi  [which.max(Yield)]

c(mmax,xmax)
 
abline(v=mmax,lty=3)
legend("topright",lwd=c(1,2),legend=c("Yield","Biomass"))


r<-2.61
K <- 1.34e8
Q<-3.8e-5
Estar <- 0.5*r/Q
Estar

writelabel("B")

par(mar=c(5.1,4.1,4.1,2.1))
