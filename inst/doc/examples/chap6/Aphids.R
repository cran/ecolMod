###############################################
## Soetaert and Herman (2008)                ##
## A practical guide to ecological modelling ##
## Chapter 6.                                ##
## Model solution - numerical methods        ##
## Chapter 6.6.4                             ##
## Aphids on a row of plants                 ##
###############################################

#----------------------#
# the model equations: #
#----------------------#
  
model <-function(t,APHIDS,parameters)
 {
    deltax     <- c (0.5,rep(1,numboxes-1),0.5)    
    Flux       <- -D*diff(c(0,APHIDS,0))/deltax
    dAPHIDS    <- -diff(Flux)/delx  + APHIDS*r

    # the output
      list(dAPHIDS )
  }  # end of model

#-----------------------#
# the model parameters: #
#-----------------------#

D         <- 0.3    # m2/day  diffusion rate
r         <- 0.01   # /day    net growth rate
delx      <- 1      # m       thickness of boxes
numboxes  <- 60 


# distance of boxes on plant, m
Distance        <- seq(from=0.5,by=delx,length.out=numboxes)  # sequence, 1 m intervals

#--------------------------#
# Initial conditions:      #
#--------------------------#

APHIDS          <- rep(0,times=numboxes)       # ind/m2   The aphid density
APHIDS[30:31]   <- 1
state           <- c(APHIDS=APHIDS)            # the state variables are initialised
                  
#----------------------#
# RUNNING the model:   #
#----------------------#

require(deSolve)
times     <-seq(0,200,by=1)   # output wanted at these time intervals           
out       <- ode.band(state,times,model,parms=0,nspec=1)  

# the data in 'out' consist of: 1st col times, 2-41: the density
# select the density data
DENSITY   <- out[,2:(numboxes  +1)]

#------------------------#
# PLOTTING model output: #
#------------------------#

# 1. temporal-spatial plot of the densities
windows()
par(oma=c(0,0,3,0))   # set outer margin size (oma)
color= topo.colors

filled.contour(x=times,y=Distance,DENSITY,color= color,
               xlab="time, days", ylab= "Distance on plant, m",main="Density")
mtext(outer=TRUE,side=3,"Aphid model",cex=1.5)  # margin text

# 2. plot initial, intermediate and final densities, density versus time..
windows()
par(mfrow=c(2,2),oma=c(0,0,3,0))   # multiple figures on a row (2 rows, 2 cols), change margin size (oma)

plot(Distance,DENSITY[1,]  ,type="l",lwd=2,xlab="Distance, m",ylab="Density", main="initial condition")
plot(Distance,DENSITY[100,],type="l",lwd=2,xlab="Distance, m",ylab="Density", main="100 days")
plot(Distance,DENSITY[200,],type="l",lwd=2,xlab="Distance, m",ylab="Density", main="200 days")

meanAphid <- rowMeans(out[,2:ncol(out)])
plot(times,meanAphid  ,type="l",lwd=2,xlab="time, days",ylab="/m2", main="Density versus time") 

mtext(outer=TRUE,side=3,"Aphid model",cex=1.5)


