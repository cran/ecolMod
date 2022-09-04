###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.7.3                             ##
## Rain of organic matter in the ocean       ##
###############################################

# load package with the integration routine:

require(deSolve)

#----------------------#
# the model equations: #
#----------------------#

model<-function(t,Conc,parameters)
 {
     Flux       <- c(flux, v*Conc)
     dConc      <- -diff(Flux)/delx  - k*Conc
     list(dConc)     # result

  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

v     <- 50.0   # m/day  sinking rate
k     <- 0.2    # /day   decay rate
flux  <- 100.0  # mmolC/m2/d  flux at upper boundary

delx      <- 10   # m      thickness of boxes
numboxes  <- 40 

# depth to air-sea interface of centre of boxes, m
Depth  <- seq(from=5,by=delx,length.out=numboxes)  # sequence, 1 m intervals

#--------------------------#
# Initial conditions:      #
#--------------------------#

state       <- rep(0,times=numboxes)
                  
#----------------------#
# RUNNING the model:   #
#----------------------#

times     <-seq(0,100,by=1)   # output wanted at these time intervals                         
out       <-ode(state,times,model,parms=0)  

# the data in 'out' consist of: 1st col times, 2-41: the concentrations
# Select the concentration data
CONC      <- out[,2:(numboxes+1)]

#------------------------#
# PLOTTING model output: #
#------------------------#

# 1. temporal-spatial plot of the concentrations
par(oma=c(0,0,3,0))   # set margin size (oma)
col <-   topo.colors
#col <-   greycol
filled.contour(x=times,y=Depth,CONC,color= col,ylim=c(395,5),
               xlab="time, days", ylab= "Depth, m",main="Concentration, mmolC/m3")
mtext(outer=TRUE,side=3,"Sinking model",cex=1.5)

# 2. Compare numeric approximation of steady-state with analytical solution
# select numerical solution (last row of CONC)
Numeric  <- CONC[nrow(CONC),]    

coords<- c(0.4,0.8,0.35,0.9)
par(fig=coords, new=TRUE)

# The analytical solution
Analytic <-  flux/v*exp(-k/v*Depth)

plot(Numeric,Depth,ylim=c(400,0),pch=21,main="steady-state",
     xlab="conc, mmolC/m3",ylab="Depth, m",cex=1.5)
lines(Analytic,Depth,lwd=2)
legend("bottomright",c("numeric","analytic"),
pch=c(21,NA),lty=c(NA,1),lwd=c(NA,2))

