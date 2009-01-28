###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 7.                                ##
## Stability and steady-state                ##
## Chapter 7.8.3                             ##
## The Lorenz equations-chaos                ##
###############################################

# package with the integration and 3-d plotting routine:

require(deSolve)
require(scatterplot3d)

#----------------------#
# the model equations: #
#----------------------#

Lorenz<-function(t,state,parameters)
  {
  with(as.list(c(state)),{ 

    ## the rates of change of state variables
    dx     <- -8/3*x+y*z
    dy     <- -10*(y-z)
    dz     <- -x*y+28*y-z

    ## the output
    list(c(dx,dy,dz))            })
 }  # end of model

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(x=1,
              y=1,
              z=1)

#----------------------#
# RUNNING the model:   #
#----------------------#

times <-seq(0,100,0.001)
out   <-as.data.frame(vode(state,times,Lorenz,0))

#------------------------#
# PLOTTING model output: #
#------------------------#

#windows()       
#par(mfrow=c(2,2),oma=c(1,2,2,2))  
#plot (out$y,out$x,type="l")
#plot (out$y,out$z,type="l")
#plot (out$z,out$x,type="l")
#plot (out$time,out$x,type="l")
#mtext(outer=TRUE,side=3,"Lorenz butterfly",line=-1)

windows()
scatterplot3d(out$x,out$y,out$z,type="l",main="Lorenz butterfly",ylab="",
              grid=FALSE,box=FALSE)

