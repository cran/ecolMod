###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.8.1                             ##
## Multiple stable states:                   ##
## the spruce budworm model                  ##
###############################################

# load package with the steady-state routines:
require(rootSolve)
require(shape)

windows(9,5)
par(mfrow=c(1,2))
ri   <- 0.05
K    <- 10    
beta <- 0.1   
ks   <- 1     

rate <- function(B, ri=0.05) 
        ri*B*(1-B/K)-beta*B^2/(B^2+ks^2)

Bsq  <- seq(0,10,length=500)
rsq  <- seq(0.01,0.07,by=0.02)
mat  <- outer(X=Bsq,Y=rsq,function(X,Y) rate(B=X,ri=Y))

matplot(Bsq,mat,xlab="B",ylab="dB/dt",
        type="l",lty=1:4,col=1) 
abline(h=0,lty=2)               
legend("bottomleft",legend=rsq,title="ri",col=1,
       lty=1:4)
writelabel("A")       

curve(rate(x,0.05),xlab="B", ylab="dB/dt",main="ri=0.05",
      from=0,to=10)
abline(h=0)


require(rootSolve)

equilibrium <- function(ri)
{
 Eq     <- uniroot.all(f=rate,interval=c(0,10),ri=ri)
 eqtype <- vector(length=length(Eq))
 for (i in 1:length(Eq))
  {  
     jac       <- gradient(f =rate,x=Eq[i],ri=ri)
     eqtype[i] <- sign(jac)+2
  }     
  return(list(x=Eq, type=eqtype))
}


eq   <- equilibrium (ri=0.05)
points(x=eq$x,y=rep(0,length(eq$x)),pch=21,cex=2,
       bg=c("grey","black","white")[eq$type] )
writelabel("B")       

windows()
rseq <- seq(0.01,0.07,by=0.0001)

plot(0,xlim=range(rseq),ylim=c(0,10),type="n",
xlab="ri",ylab="B*",main="spruce budworm model")

for (ri in rseq) {
eq  <- equilibrium(ri) 

points(rep(ri,length(eq$x)),eq$x,pch=22,
col=c("darkgrey","black","lightgrey")[eq$type],
bg =c("darkgrey","black","lightgrey")[eq$type]) 
}

equi <-uniroot.all(f=rate,interval=c(0,10),r=0.05)
arrows(0.05,10         ,0.05,equi[4]+0.2,length=0.1 )
arrows(0.05,equi[3]+0.2,0.05,equi[4]-0.2,length=0.1 )
arrows(0.05,equi[3]-0.2,0.05,equi[2]+0.2,length=0.1 )
arrows(0.05,equi[1]+0.2,0.05,equi[2]-0.2,length=0.1 )  

equi <-uniroot.all(f=rate,interval=c(0,10),r=0.038) 
arrows(0.038,10         ,0.038,equi[2]+0.1,length=0.1 )
arrows(0.038,equi[1]+0.1,0.038,equi[2]-0.1,length=0.1 )

equi <-uniroot.all(f=rate,interval=c(0,10),r=0.07) 
arrows(0.07,10         ,0.07,equi[2]+0.2,length=0.1 )
arrows(0.07,equi[1]+0.2,0.07,equi[2]-0.2,length=0.1 )
