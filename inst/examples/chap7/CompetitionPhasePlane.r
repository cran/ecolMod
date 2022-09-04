###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.8.5                             ##
## Fate of marine zooplankton in an estuary  ##
## equilibrium condition                     ##
###############################################

# load package with the steady-state and integration routines:
require(rootSolve)
require(deSolve)

windows()
r1    <- 3              # parameters
r2    <- 2
K1    <- 1.5
K2    <- 2
alf12 <-1
alf21 <-2

Lotka<-function(t,N,pars)

{
 dN1 <- r1*N[1]*(1-(N[1]+alf12* N[2])/K1)
 dN2 <- r2*N[2]*(1-(N[2]+alf21* N[1])/K2)

 list(c(dN1  , dN2  ))     # the rate of change
}

Ax  <- c(0,K2/alf21)
Ay  <- K2 - alf21* Ax
By  <- c(0,K1/alf12)
Bx  <- K1 - alf12* By
xlim   <- range(c(Ax, Bx))
ylim   <- range(c(Ay, By))

plot  (x=Ax,y=Ay, type="l", lwd=3,   # 1st isocline
     main="Competition phase-plane",
       xlab="N1",ylab="N2",xlim=xlim,ylim=ylim)
lines (Bx,By,lwd=3,lty=2)            # 2nd isocline

tarrows <- function(out,ds,...)
{
 # select point at predefined distance from point (p1,p2)
 p   <- unlist(out[nrow(out),2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i1<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)

 p   <- unlist(out[1,2:3])
 dd  <- (out[,2]-p[1])^2+ (out[,3]-p[2])^2
 dd2 <-  c(dd[-1],dd[1])
 i2<-which(dd<ds&dd2>ds | dd>ds&dd2<ds)[1]
 ii <- c(i1,i2)
# ii <- which(dd>ds)
# iseq <- seq(1,length (ii),15)
 for (i in ii ) arrows(out[i,2],out[i,3],out[i+1,2],out[i+1,3],length=0.1,lwd=1,...)

}

trajectory <- function(N1,N2)
{
times  <-seq(0,30,0.1)
state  <-c(N1 = N1, N2 = N2)
out    <-as.data.frame(ode(state,times,Lotka,0))

lines (out$N1,out$N2,type="l")
tarrows(out,ds)
}

ds<- 0.1

trajectory (0.05,0.3)
trajectory (0.11,0.3)
trajectory (1.5,1.8)
trajectory (1.0,2.0)

# 4 equilibrium points
X  <- c(0,0 ,K1,(K1-alf12*K2)/(1-alf12*alf21))
Y  <- c(0,K2,0 ,(K2-alf21*K1)/(1-alf12*alf21))

# Jacobian matrix, and eigenvalues
ei    <- matrix(nrow=4,ncol=2)

for ( i in 1:4)   # for each root
{
 # jacobian matrix
 N1 = X[i]
 N2 = Y[i]
 Jacob <- jacobian.full(y=c(N1,N2),func=Lotka)

 # eigenvalues
 ei[i,]   <- eigen(Jacob)$values

 # white:unstable node, black:stable node, grey:saddle
 if (sign(ei[i,1])>0 & sign(ei[i,2])>=0) col <- "white"
 if (sign(ei[i,1])<0 & sign(ei[i,2])<=0) col <- "black"
 if (sign(ei[i,1])* sign(ei[i,2])   <0 ) col <- "grey"

# equilibrium point plotting
 points(N1,N2,pch=22,cex=2.0,bg=col,col="black")
}
cbind(N1=X,N2=Y,ei)

eig      <- eigen(Jacob)
vv       <- eig$vector[eig$values<0]



# the reverse model, output -rate of change
revmod <- function(t,N,p) list(-1*unlist(Lotka(t,N,p)))
times  <-seq(0,1.9,0.05)

# first direction
state <- c(N1,N2) + 0.01*vv

out   <-as.data.frame(ode(state,times,revmod,0))
lines(out[,2],out[,3],lty=2)
tarrows(out,ds,code=1)

# second direction
state <- c(N1,N2) - 0.01*vv
times  <-seq(0,10,0.05)
out   <-as.data.frame(ode(state,times,revmod,0))
lines(out[,2],out[,3],lty=2)
tarrows(out,ds,code=1)

# trajectories out of the equilibrium point
ww       <- eig$vector[eig$values>0]
trajectory (N1+0.05*ww[1],N2+0.05*ww[2])
trajectory (N1-0.05*ww[1],N2-0.05*ww[2])

legend("right",legend=c("isocline N1","isocline N2","trajectory","separatrice",
"saddle point","stable equilibrium","unstable equilibrium"),lty=c(2,1,1,2,NA,NA,NA),
lwd=c(2,2,1,1,NA,NA,NA),pch= c(NA,NA,NA,NA,22,22,22),pt.bg=c(NA,NA,NA,NA,"grey","black","white"))

