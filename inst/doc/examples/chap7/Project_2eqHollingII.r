###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.9.4. Project                    ##
## Predator-prey system with type II         ##
## functional response                       ##
###############################################

library(rootSolve)
library(deSolve)

#=======================================================================
# Stability properties of model with density dependent growth of prey
# and holling type II functional response of grazing
#=======================================================================

# parameters
r     <- 3      # prey rate of increase        
K     <- 10     # prey carrying capacity
ks    <- 2      # half-saturation constant
e     <- 0.5    # growth efficiency of predator
m     <- 0.2    # predator mortality rate
In    <- 0.6    # ingestion rate predator

#========================================
# rate of change 
#========================================
model<-function(t,State,pars)
{
 with (as.list(State),{
 dN <- r*Prey*(1-Prey/K)-In*Prey/(Prey+ks)*Predator
 dP <- e*In*Prey/(Prey+ks)*Predator - m*Predator

 list(c(dN  , dP  ))     # the rate of change
 })
}

#========================================
# isoclines and roots
#========================================
isocline <- function()
{

# prey isocline is a curve
 curve  (r/In*(1-x/K)*(x+ks), from=0, to= K,lwd=3,   # 1st isocline
        xlab="Prey",ylab="Predator",xlim=xlim,ylim=ylim)

# predator isocline is a vertical line
 Neq    <- m*ks/(In*e-m)        
 if (Neq < 0 | is.infinite(Neq)) {warnings("Model has no solution");return}
 abline (v=Neq,lwd=3,lty=2)                # 2nd isocline
# legend("topright",legend=paste("r=",In))

# first equilibrium point : (0,0)
# then intersection between predator isocline and X-axis, 
# and intersection between predator and prey isocline
 roots <- matrix(ncol=2,byrow=TRUE,data=
          c(0  ,0,
            K  ,0,
            Neq, r/In*(1-Neq/K)*(Neq+ks)))
 colnames(roots) <- c("Prey","Predator")

 root<-stability(roots)
 return(root)
 
}
# end isocline

#========================================
# the stability properties of the equilibrium points
#========================================

stability <- function (roots)
{
 small    <- 1e-8
 # Jacobian matrix, and eigenvalues
 eig    <- NULL
 for (i in 1:nrow(roots))
 {
  equi <- roots[i,]
  Jacob <- jacobian.full(y=equi,func=model)

 # eigenvalues
  ei   <- eigen(Jacob)$values
  eig <- rbind(eig,ei)

 # white:unstable node, black:stable node, grey:saddle
  ifelse (is.complex(ei), pch<-21, pch<-22)
  ei <- Re(ei)
  if (sign(ei[1])>0 & sign(ei[2])>=0) pcol <- "white"
  if (sign(ei[1])<0 & sign(ei[2])<=0) pcol <- "black"
  if (sign(ei[1])   * sign(ei[2])<0 ) pcol <- "grey"
 # equilibrium point plotting

  points (equi[1],equi[2],cex=2,pch=pch,bg=pcol,col="black")

 }
 return(list(equilibrium=roots,eigenvalues=eig) )
}   # end stability


#========================================
# model trajectories
#========================================

trajectory <- function(N,P)
{
times  <-seq(0,100,0.1)
state  <-c(Prey = N, Predator = P)
out    <-as.data.frame(ode(state,times,model,0))

lines (out$Prey,out$Predator,type="l")
arrows(out[10,2],out[10,3],out[11,2],out[11,3],length=0.1,lwd=1)

}

#========================================
# model applications
#========================================
windows()
par(oma=c(0,0,1,0),mfrow=c(2,2))

In     <- 0.55    # ingestion rate predator
xlim   <- c(0,K)
ylim   <- c(0,25)
isocline( ) 
trajectory (6,10)
trajectory (10,15)
trajectory (4,5)
trajectory (2,20)
title("stable equilibrium")

In     <- 0.45    # ingestion rate predator
ylim   <- c(0,30)
xlim   <- c(0,20)
isocline( ) 
trajectory (6,5)
trajectory (10,25)
trajectory (20,20)

title("stable equilibrium")

In     <- 0.8    # ingestion rate predator
ylim   <- c(0,20)
xlim   <- c(0,K)
isocline( ) 
trajectory (6,8)
trajectory (10,5)
trajectory (1,1)
trajectory (2.1,12.1)
title("stable limit cycle")

mtext(3,text="Predator-Prey, Monod grazing",line=-1,outer=TRUE,cex=1.2)


plot(0,type="n",axes=FALSE,xlab="",ylab="")
legend("center",pch=c(22,22,21,NA,NA),pt.bg=c("black","grey","white",NA,NA),
      lty=c(NA,NA,NA,1,2),lwd=c(NA,NA,NA,2,2),pt.cex=2,title="equilibrium",
      legend=c("stable","unstable","neutral","constant prey","constant predator"))

