###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.6.1                             ##
## The enzymatic reaction model              ##
###############################################

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,state,parameters){
 with(as.list(c(state,parameters)),{  # unpack the state variables, parameters

    dD <- -k1*E*D + k2*I
    dI <-  k1*E*D - k2*I - k3*I*F
    dE <- -k1*E*D + k2*I + k3*I*F
    dF <-                - k3*I*F
    dG <-                  k3*I*F
    list(c(dD,dI,dE,dF,dG))          # the output, packed as a list
    })
}


#-----------------------#
# the model parameters: #
#-----------------------#
                                     
parameters<-c(k1=0.01/24,            # parameter values
              k2=0.1/24,
              k3=0.1/24)

#-------------------------#
# the initial conditions: #
#-------------------------#
 
state     <-c(D=100,                # state variable initial conditions
              I=10,
              E=1,
              F=1,
              G=0)

#----------------------#
# RUNNING the model:   #
#----------------------#

times     <-seq(0,300,0.5)         # time: from 0 to 1000 hours, steps of 0.5 hours

# integrate the model 
require(deSolve)                   # package with the integration routine:
out <-as.data.frame(ode(state,times,model,parameters))

#------------------------#
# PLOTTING model output: #
#------------------------#
#windows()
par(mfrow=c(2,2), oma=c(0,0,3,0))   # set number of plots (mfrow) and margin size (oma)

plot (times,out$D,type="l",main="[D]",xlab="time, hours",ylab="mol/m3",lwd=2)
plot (times,out$F,type="l",main="[F]",xlab="time, hours",ylab="mol/m3",lwd=2)
plot (times,out$E,type="l",main="[E]",xlab="time, hours",ylab="mol/m3",lwd=2)
plot (times,out$I,type="l",main="[I]",xlab="time, hours",ylab="mol/m3",lwd=2)
mtext(outer=TRUE,side=3,"enzymatic reaction",cex=1.5)
#plot (times,out$sum,type="l",col="red")
